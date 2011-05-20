/**  @file BkgndTupleSelectTool.cxx
    @brief implementation of class BkgndTupleSelectTool
    
  $Header$  
*/

#include "IBkgndTupleSelectTool.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IDataProviderSvc.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/MCEvent.h"

#include "facilities/Util.h"

#include "BackgroundManager.h"
#include "XmlFetchEvents.h"
#include "IInterleaveMap.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLeaf.h"

#include "CLHEP/Random/RandFlat.h"

// access to the tuple
#include "ntupleWriterSvc/INTupleWriterSvc.h"
// access to the CEL
#include "RootIo/IRootIoSvc.h"

#include <stdexcept>
#include <cassert>
#include <sstream>
#include <cmath>
#include <vector>
#include <list>
#include <ctime>

/** @class BackgroundSelection
    @brief manage the selection of background events to merge with signal events
    @author Dan Flath

It uses the name of a tuple variable as a key to look up trigger and downlink rates of a background
source depending on the value of that variable. See the corresponding member functions. Given such
a value, it expects to find a source of such events generated with that value which can be inserted 
into the Gleam output tuple.

*/
class BkgndTupleSelectTool : public AlgTool, virtual public IBkgndTupleSelectTool
{
public:

    // Standard Gaudi Tool constructor
    BkgndTupleSelectTool(const std::string& type, const std::string& name, const IInterface* parent);

    // After building it, destroy it
    ~BkgndTupleSelectTool();

    /// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief Finalize method for the tool
    StatusCode finalize();

    ///! The current value of the quantity that we are selecting on
    double value()const;

    /** @brief select an event and copy the contents to the output tree
    */
    void selectEvent();

    virtual void copyEventInfo();

    /**@brief the current trigger rate for the current value
    */
    double triggerRate()const{return m_triggerRate;}

    /**@brief the downlink rate for the current value
    */
    double downlinkRate()const{return m_downlinkRate;}

    const std::string& sourceName() const;
    //const std::string& name() const {return this->name();}

private:
    /** @brief disable branches in the tree, from the disableList
    */
    void disableBranches(TTree*);

    /**@brief Set the addresses of active leaves, so that a GetEvent will perform a copy
    */
    void setLeafPointers();

    /**@brief Retrieve the correct tree for current value of the variable(s) so copyEvent comes from correct sample bin
    */
    void setCurrentTree(double val);

    /// access the ntupleWriter service to write out to ROOT ntuples
    INTupleWriterSvc *   m_rootTupleSvc;

    /// access the RootIoSvc to get the CompositeEventList ptr
    IRootIoSvc *         m_rootIoSvc;

    TLeaf*               m_varLeaf;         ///< corresponding TLeaf for access to current value
    TLeaf*               m_runLeaf;
    TLeaf*               m_eventLeaf;
    TLeaf*               m_timeLeaf;
    TLeaf*               m_liveLeaf;
    TLeaf*               m_mcidLeaf;

    TTree*               m_outputTree; ///< the tree to copy data to
    TTree*               m_inputTree; ///< the tree to copy data from

    unsigned int         m_eventOffset;

    std::string          m_tupVarName;

    IFetchEvents*        m_fetch; ///< abstract guy that processes the xml file
    double               m_downlinkRate;
    double               m_triggerRate;

    // For dealing with the interleave <-> output ntuple mapping
    IInterleaveMap*      m_interleaveMap;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*    m_dataSvc;

    std::list<double>    m_binList;

    StringProperty       m_treeName;    ///< name of the tree to process
    StringArrayProperty  m_disableList;
};

static ToolFactory<BkgndTupleSelectTool> s_factory;
const IToolFactory& BkgndTupleSelectToolFactory = s_factory;

//------------------------------------------------------------------------
BkgndTupleSelectTool::BkgndTupleSelectTool(const std::string& type, 
                                 const std::string& name, 
                                 const IInterface* parent) :
                                 AlgTool(type, name, parent)
                                 , m_rootTupleSvc(0)
                                 , m_rootIoSvc(0)
                                 , m_varLeaf(0)
                                 , m_runLeaf(0)
                                 , m_eventLeaf(0)
                                 , m_outputTree(0)
                                 , m_inputTree(0)
                                 , m_eventOffset(0)
                                 , m_fetch(0)
                                 , m_downlinkRate(0)
                                 , m_triggerRate(0)
{
    //Declare the additional interface
    declareInterface<IBkgndTupleSelectTool>(this);

    // declare properties with setProperties calls
    declareProperty("TreeName",       m_treeName="MeritTuple");
    declareProperty("DisableList",    m_disableList);

    std::string varName = this->name();
    int         dotPos  = varName.find(".");
    m_tupVarName = varName.substr(dotPos+1,varName.size());
    m_binList.clear();
}
//------------------------------------------------------------------------
BkgndTupleSelectTool::~BkgndTupleSelectTool()
{
    delete m_fetch;
    delete m_inputTree;
}

StatusCode BkgndTupleSelectTool::initialize()
{
    StatusCode sc   = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // Set the properties
    setProperties();

    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        log << MSG::ERROR << "Service [EventDataSvc] not found" << endreq;
        return sc;
    }

    // Get the interleave map tool
    IAlgTool* toolPtr = 0;
    if ( (sc = toolSvc()->retrieveTool("InterleaveMapTool", toolPtr)).isFailure())
    {
        log << MSG::ERROR << "failed to retrieve InterleaveMapTool for " + name() << endreq;
        return sc;
    }
    else m_interleaveMap = dynamic_cast<IInterleaveMap*>(toolPtr);

    // Retrieve tuple info an set up
    // get a pointer to RootTupleSvc
    sc = service("RootTupleSvc", m_rootTupleSvc);       
    if( sc.isFailure() ) 
    {
        log << MSG::ERROR << "failed to get the RootTupleSvc" << endreq;
        return sc;
    }

    // Use this to recover a pointer to the output tree
    void * ptr= 0;
    m_rootTupleSvc->getItem(m_treeName.value().c_str(),"", ptr);
    if( ptr==0)
    {
        log << MSG::ERROR << "Could not find the tuple: " +  m_treeName.value() << endreq;
        return StatusCode::FAILURE;
    }
    m_outputTree = static_cast<TTree*>(ptr);

    sc = service("RootIoSvc", m_rootIoSvc, true);       
    if( sc.isFailure() ) 
    {
        log << MSG::ERROR << "failed to get the RootIoSvc" << endreq;
        return sc;
    }
  
////    BackgroundManager::instance()->getCelManager()->addComponent(m_treeName.value(), m_outputTree);

    std::string  filePath;

    //    StringProperty filePathProperty("FilePath", filePath);

    // Retrieve the file path from the parent algorithm (which can be tricky)
    // The try-catch is for the services (e.g. FluxSvc) which instantiate all
    // tools and for which we fail to find the desired property
    try
    {
        if (IProperty* parentProp = dynamic_cast<IProperty*>(const_cast<IInterface*>(parent())))
        {
            SimplePropertyRef<std::string> filePathProperty("FilePath", filePath);
            StatusCode sc = parentProp->getProperty(&filePathProperty);
            sc = parentProp->getProperty(&m_disableList);
        }
    }
    catch(...) 
    {
        int j = 0;
    }

    if( !filePath.empty())
    {
        facilities::Util::expandEnvVar(&filePath);
        log << MSG::INFO << "Using xml file path " << filePath 
            << " for BkgndTupleSelectTool " + name() + "." << endreq;

        // If passed an XML file, set up the XmlFetchEvents object
        std::string xmlFileName = filePath+"/" + m_tupVarName + ".xml";

        m_fetch = new XmlFetchEvents(xmlFileName, m_tupVarName);

        // find leaf corresponding to varname in the output tree
        if( (m_varLeaf = m_outputTree->GetLeaf(m_tupVarName.c_str())) == 0 ) 
        {
            log << MSG::ERROR << "Could not find leaf pointer for " + name() << endreq;
            return StatusCode::FAILURE;
        }
        if( (m_runLeaf = m_outputTree->GetLeaf("EvtRun")) == 0)
        {
            log << MSG::ERROR << "EvtRun leaf not found in the tuple" << endreq;
            return StatusCode::FAILURE;
        }
        if( (m_eventLeaf = m_outputTree->GetLeaf("EvtEventId")) == 0)
        {
            log << MSG::ERROR << "EvtEventId leaf not found in the tuple" << endreq;
            return StatusCode::FAILURE;
        }
        if( (m_timeLeaf = m_outputTree->GetLeaf("EvtElapsedTime")) == 0)
        {
            log << MSG::ERROR << "EvtElapsedTime leaf not found in the tuple" << endreq;
            return StatusCode::FAILURE;
        }
        if( (m_liveLeaf = m_outputTree->GetLeaf("EvtLiveTime")) == 0)
        {
            log << MSG::ERROR << "EvtLiveTime leaf not found in the tuple" << endreq;
            return StatusCode::FAILURE;
        }
        if( (m_mcidLeaf = m_outputTree->GetLeaf("McSourceId")) == 0)
        {
            log << MSG::ERROR << "McSourceId leaf not found in the tuple" << endreq;
            return StatusCode::FAILURE;
        }

        // This is unnecessary and is being eliminated
        //// prime the pump with the mean of the limits for initial downlink, trigger rates
        //double mean( (m_fetch->minValFullRange() + m_fetch->maxValFullRange())/2.);

        //setCurrentTree(mean);
    }

    return sc;
}

StatusCode BkgndTupleSelectTool::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;
    
    return status;
}

const std::string& BkgndTupleSelectTool::sourceName() const 
{
    return m_fetch->name();
}

//------------------------------------------------------------------------
double BkgndTupleSelectTool::value()const
{
    double x(m_varLeaf->GetValue());
    return x;
}
//------------------------------------------------------------------------
void BkgndTupleSelectTool::selectEvent()
{
    MsgStream log(msgSvc(), name());
    try {
    double x(value());

#if 0 // have to do this since my specification was not understood
    x = fabs(x); 
#endif

    if( !m_fetch->isValid(x) )
    {
        std::stringstream msg;
        msg << "BkgndTupleSelectTool::selectEvent: called with " << name() 
            <<" = "<< x << " is not in range " 
            <<  m_fetch->minValFullRange() << ", to " 
            <<  m_fetch->maxValFullRange() << std::endl;

        throw std::runtime_error(msg.str());
    }

    // make sure we have the right tree selected for new value    
    // if still valid, do not change
    if( !m_fetch->isCurrent(x) )
    {
        setCurrentTree(x);   
    }

    // grab the next event, cycling if needed
    int code = m_inputTree->GetEntry(m_eventOffset++);

    if( code <= 0)
    { 
        m_eventOffset = 0; 
        time_t rawtime;
        time(&rawtime);
        log << MSG::DEBUG << ctime(&rawtime) << "Failed to GetEntry, "
            << "trying event 0" << endreq;
        code = m_inputTree->GetEntry(m_eventOffset++);
        if( code <= 0 )
        {
            throw std::runtime_error("BkgndTupleSelectTool::selectEvent -- could not read file");
        }
    }
    } catch(...) {
        log << MSG::WARNING << "exception thrown" << endreq;
        throw;
    }
}
//------------------------------------------------------------------------
void BkgndTupleSelectTool::disableBranches(TTree* t) 
{
    const std::vector<std::string >& disableList = m_disableList;
    for( std::vector<std::string >::const_iterator it= disableList.begin(); it != disableList.end();  ++it) 
    {
        t->SetBranchStatus((*it).c_str(), false);
    }

}
//------------------------------------------------------------------------
void BkgndTupleSelectTool::setCurrentTree(double x) 
{
    MsgStream log(msgSvc(), name());

   try {
    // replace the TChain
    delete m_inputTree; 

    std::string treeName = m_treeName;
    m_inputTree = new TChain(treeName.data());

    // this is necessary due to the design of ROOT :-(
    TDirectory *saveDir = gDirectory;
    
    int stat = m_fetch->getFiles(x, dynamic_cast<TChain*>(m_inputTree), log.level()<=MSG::DEBUG);
    
    if( stat!=0 )
    {
        std::stringstream msg;
        msg << "BackgaroundSelection::setCurrentTree: called with " << name() 
            << " = " << x << ", getFiles returned error code " << stat;

        throw std::runtime_error(msg.str());
    }

    // this is necessary due to the design of ROOT :-(
    //saveDir->cd(); HMK Move this to the end of the method

    // start at a random location in the tree:
    double length (m_inputTree->GetEntries());
    if( length==0 ) 
    {
        throw std::runtime_error("BkgndTupleSelectTool::setCurrentTree: no events in the tree");
    }
    
    m_eventOffset = (unsigned int)(CLHEP::RandFlat::shoot()*(length - 1));

    // point tree to buffer for copying events:
    setLeafPointers();
    m_triggerRate  = m_fetch->getAttributeValue("triggerRate", x);
    m_downlinkRate = m_fetch->getAttributeValue("downlinkRate", x);
    log << MSG::DEBUG << "triggerRate: " << m_triggerRate 
        << " downlinkRate: " << m_downlinkRate << endreq;

    int ret = m_inputTree->GetEntry(0);
    if( ret<0 ) throw std::runtime_error("BkgndTupleSelectTool::setCurrentTree: could not read");

    double minBin = m_fetch->minVal();

    std::list<double>::iterator binIter = std::find(m_binList.begin(),m_binList.end(),minBin);
    if (binIter != m_binList.end())
    {
        m_binList.push_back(minBin);

////        BackgroundManager::instance()->getCelManager()->addComponent(m_treeName.value(), m_inputTree);
    }

    // this is necessary due to the design of ROOT :-(
    saveDir->cd();
    } catch( const std::exception& e) {
        log << MSG::WARNING << e.what() << endreq;
        throw(e);
   } catch(...) {
      log << MSG::WARNING << "exception thrown" << endreq;
      throw;
   }

    return;
}
//------------------------------------------------------------------------
void BkgndTupleSelectTool::setLeafPointers()
{
    // first disable the branches
    disableBranches(m_inputTree);

    // iteration over active leaves -- copied from Ttree::Show()
    TObjArray *leaves  = m_inputTree->GetListOfLeaves();
    if( leaves==0 ) throw std::runtime_error("BkgndTupleSelectTool::setLeafPointers: invalid tree, no leaves found");

    Int_t nleaves = leaves->GetEntries();
    if (nleaves==0) throw std::runtime_error("BkgndTupleSelectTool::setLeafPointers: invalid tree, no leaves found");
    for (Int_t i=0;i<nleaves;i++) 
    {
        TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(i);
        TBranch *branch = leaf->GetBranch();
        if (branch->TestBit(kDoNotProcess)) continue;
        Int_t len = leaf->GetLen();
        if (len <= 0) continue;

        // here copy the value pointer from the output tuple
        TLeaf* otherleaf = m_outputTree->GetLeaf(leaf->GetName());
        if( otherleaf==0)
        {
            std::stringstream msg;
            msg << "BkgndTupleSelectTool::setLeafPointers: did not find the leaf " << leaf->GetName();
            std::cerr << msg.str() << std::endl;
            throw std::runtime_error( msg.str());
        }
        else
        {
            m_inputTree->SetBranchAddress(leaf->GetName(), otherleaf->GetValuePointer());
        }
    }
}
//------------------------------------------------------------------------
namespace {
    //! Utility to set a scalar value in a ROOT tree. 
    //! Assume that the type is known!
    template <typename  Type> 
        void setLeafValue(TLeaf* leaf, Type newvalue)
    {
        const char* cLeafType = leaf->GetTypeName();
        std::string leafType(cLeafType);
        void*       leafPtr   = leaf->GetValuePointer();

        if (leafType == "UInt_t")
        {
            *(reinterpret_cast<UInt_t*>(leafPtr)) = newvalue;
        }
        else if (leafType == "Float_t")
        {
            *(reinterpret_cast<Float_t*>(leafPtr)) = newvalue;
        }
        else if (leafType == "Double_t")
        {
            *(reinterpret_cast<Double_t*>(leafPtr)) = newvalue;
        }
        else
        {
            *(reinterpret_cast<Type*>(leafPtr)) = newvalue;
        }
    }

}

void BkgndTupleSelectTool::copyEventInfo()
{
    SmartDataPtr<Event::EventHeader> header(m_dataSvc, EventModel::EventHeader);

    // these types *must* correspond with those in EvtValsTool, which this code replaces for interleaved events
    float  EvtRun         = header->run();
    float  EvtEventId     = header->event();
    double EvtElapsedTime = header->time();
    double EvtLiveTime    = header->livetime();

    // load the map entries to relate the interleaved row back to its original source
    m_interleaveMap->getInterleaveMap()->run      =  header->run();
    m_interleaveMap->getInterleaveMap()->curEvent =  header->event();
    m_interleaveMap->getInterleaveMap()->irun     =  static_cast<int>(m_runLeaf->GetValue());
    m_interleaveMap->getInterleaveMap()->ievent   =  static_cast<int>( m_eventLeaf->GetValue());
    m_interleaveMap->getInterleaveMap()->value    =  value();
    strncpy(m_interleaveMap->getInterleaveMap()->type, m_tupVarName.c_str(), 
        sizeof(m_interleaveMap->getInterleaveMap()->type));

    // Update this event's info 
    setLeafValue(m_runLeaf,   EvtRun);
    setLeafValue(m_eventLeaf, EvtEventId);
    setLeafValue(m_timeLeaf,  EvtElapsedTime);
    setLeafValue(m_liveLeaf,  EvtLiveTime);

    // Convention is that interleaved source ids are negative of what comes from FluxSvc
    int sourceId = -static_cast<int>(m_mcidLeaf->GetValue()); 
    setLeafValue(m_mcidLeaf, sourceId);

    // Save this row
    m_rootTupleSvc->storeRowFlag(m_treeName.value(),true);
    m_interleaveMap->saveInterleaveMapRow();

    // Tell the CEL about this...
    // Recall David's note:
    // for each tree, if the associated file is not writable (this is a tree where you just get an entry), 
    // the "current" entry is obtained through a call to tree::GetReadEntry(). On the contrary,
    // if the associated file is writable (this is a tree where you just filled
    // an entry), the "current" is assumed to be (GetEntries()-1).

    // Calling CEL's fill after RootTupleSvc::saveRow which does a TTree::Fill
    // HMK March 5, 2008 Waiting until we're sure we can insert events with NULL components
    // For now we will run an after-burner to fill the interleaved events into the CEL
    //// m_rootIoSvc->getCell()->fillEntry("merit",m_outputTree);
    // If we just fill the merit entry, will the others (mc, digi, etc) be left as NULL for this event as intended?
    ////BackgroundManager::instance()->getCel()->fillEvent();
}

