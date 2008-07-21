/** @file InterleaveAlg.cxx

@brief declaration and definition of the class InterleaveAlg

$Header$

*/
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"

#include "BackgroundManager.h"
#include "IBkgndTupleSelectTool.h"

// this service manages livetime calculations
#include "Trigger/ILivetimeSvc.h"

// TDS class declarations: input data, and McParticle tree
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/MCEvent.h"

// access to the tuple
#include "ntupleWriterSvc/INTupleWriterSvc.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TLeaf.h"

// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"

#include "facilities/Util.h"

#include <cassert>
#include <map>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

/** @class InterleaveAlg
    @brief Gaudi Alg that will fill the ntuple with a background event

    When initialized, it expects to find that one or more instances of SampledBackground have declared sources, 
    with the defineSource static method 

    The execute method is expected to be invoked for every non-triggering event, on a branch for the Trigger
    algorithm sequence.

 */


class InterleaveAlg : public Algorithm {
public:
    InterleaveAlg(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 
    /// access the ntupleWriter service to write out to ROOT ntuples
    ILivetimeSvc *       m_LivetimeSvc;

    StringProperty       m_treeName;    ///< name of the tree to process
    StringProperty       m_filePath;    ///< path to files containing info for each variable
    StringProperty       m_mapName;     ///< name of the map tree
    StringArrayProperty  m_disableList;
    int                  m_count;       ///< number of processed events

    BackgroundManager*   m_bkgndManager;
};

//------------------------------------------------------------------------

static const AlgFactory<InterleaveAlg>  Factory;
const IAlgFactory& InterleaveAlgFactory = Factory;

//------------------------------------------------------------------------
//! ctor
InterleaveAlg::InterleaveAlg(const std::string& name, ISvcLocator* pSvcLocator)
              :Algorithm(name, pSvcLocator), m_count(0)
{
    // declare properties with setProperties calls
    declareProperty("FilePath",       m_filePath="");
    declareProperty("TreeName",       m_treeName="MeritTuple");
    declareProperty("MapName",        m_mapName="interleave_map");
    declareProperty("DisableList",    m_disableList);

    // initialize the disable list, which can be added to, or replaced in the JO
    std::vector<std::string> tlist;
    tlist.push_back("EvtElapsedTime");
    tlist.push_back("EvtLiveTime");
    tlist.push_back("Pt*");   // use current position, etc.
    tlist.push_back("FT1*");  // will recalculate
    m_disableList = tlist;

    // Get the pointer to the background manager
    m_bkgndManager = BackgroundManager::instance();
}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode InterleaveAlg::initialize()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    sc = service("LivetimeSvc", m_LivetimeSvc);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "failed to get the LivetimeSvc" << endreq;
        return sc;
    }

    // Make sure the file path has been set...
    if(m_filePath.value().empty())
    {
        log << MSG::INFO<< "No file specified for interleave" << endreq;
        return sc;
    }
    else
    {
        std::string filePath(m_filePath.value());
        facilities::Util::expandEnvVar(&filePath);
        log << MSG::INFO << "Using xml file path " << filePath << " for interleave." << endreq;
    }

    // Initialize the CEL manager
////    std::string celFileName = "$GLEAMDATAPATH/CELInterleaveFile.root";
////    m_bkgndManager->getCelManager()->initWrite(celFileName,"RECREATE");

    // initialize the background selection
    try 
    {
        const BackgroundManager::BkgndSrcNameVec& sourcenames   = m_bkgndManager->getBkgndSrcNameVec();
        BackgroundManager::BkgndNameToSelMap&     bkgndToSelMap = m_bkgndManager->getBkgndNameToSelMap();
        BackgroundManager::BkgndNameToSelMap&     srcToSelMap   = m_bkgndManager->getSrcNameToSelMap();

        for( BackgroundManager::BkgndSrcNameVec::const_iterator it(sourcenames.begin()); it!=sourcenames.end();++it)
        {
            const std::string& name=*it;
            log << MSG::INFO << "setting up tuple key "<<  name <<endreq;    
            
            IBkgndTupleSelectTool* background = 0;
            if ( (sc = toolSvc()->retrieveTool("BkgndTupleSelectTool", name, background, this)).isFailure())
            {
                log << MSG::ERROR << "failed to retrieve BkgndSelectTool for " + name << endreq;
                return sc;
            }

            bkgndToSelMap[name]  = background;
            std::string srcName  = background->sourceName();
            srcToSelMap[srcName] = background;
            
            log <<MSG::INFO << "\tinitial trigger, downlink rates: " << bkgndToSelMap[name]->triggerRate() 
                <<", " << bkgndToSelMap[name]->downlinkRate() 
                <<", for source name " << background->sourceName() <<  endreq;
        }
    }
    catch( const std::exception& e)
    {
        log << MSG::ERROR << e.what() << endreq;

        return StatusCode::FAILURE;
    }

    log << MSG::INFO << "initialized OK: initial downlink rate to merge is " << m_bkgndManager->downlinkRate() << " Hz"<< endreq;

    m_LivetimeSvc->setTriggerRate(m_bkgndManager->triggerRate() - m_bkgndManager->downlinkRate());

    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode InterleaveAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;

    MsgStream   log( msgSvc(), name() );

    // Keep track of the number of times we are called
    ++m_count;

    // Retrieve the source name from the MC Header
    SmartDataPtr<Event::MCEvent>  mcHeader(eventSvc(), EventModel::MC::Event );

    std::string sourceName = mcHeader->getSourceName();
    
    BackgroundManager::BkgndNameToSelMap& srcToSelMap = m_bkgndManager->getSrcNameToSelMap();
    
    if( srcToSelMap.find(sourceName) == srcToSelMap.end() )
    {
        log << MSG::ERROR << "Sourcename " << sourceName << " not found " << endreq;
        return StatusCode::FAILURE;
    }

    IBkgndTupleSelectTool* selector = srcToSelMap[sourceName];

    // ask for a tree corresponding to our current position: it will set all the tuple
    selector->selectEvent();

    // let the livetime service know about the current trigger rate,
    // and set the current accumulated live time in the header
    log << MSG::DEBUG 
        << "Event type " << selector->name()
        << " value "     << selector->value()
        << ", trigger, downlink rates: " << m_bkgndManager->triggerRate()<<", " 
        << m_bkgndManager->downlinkRate() << endreq;

    m_LivetimeSvc->setTriggerRate(m_bkgndManager->triggerRate() - m_bkgndManager->downlinkRate());
    SmartDataPtr<Event::EventHeader>   header(eventSvc(),    EventModel::EventHeader);
    header->setLivetime( m_LivetimeSvc->livetime());

    // overwrite the event info
    selector->copyEventInfo();

    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode InterleaveAlg::finalize()
{
    StatusCode  sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    log << MSG::INFO << "Inserted "<< m_count << " sampled background events." << endreq; 

////    m_bkgndManager->getCelManager()->fillFileAndTreeSet();

    return sc;
}

