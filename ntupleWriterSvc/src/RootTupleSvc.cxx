/** @file RootTupleSvc.cxx
@brief declare, implement the class RootTupleSvc
*/
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"
#include "facilities/Util.h"
#include <map>

// root includes
#include "TTree.h"
#include "TFile.h"


/** 
* @class RootTupleSvc
* @brief Special service that directly writes ROOT tuples
*
* $Header$
*/
class RootTupleSvc :  public Service, virtual public IIncidentListener,
    virtual public INTupleWriterSvc
{  

public:

    /// perform initializations for this service - required by Gaudi
    virtual StatusCode initialize ();

    /// clean up after processing all events - required by Gaudi
    virtual StatusCode finalize ();

    /// Handles incidents, implementing IIncidentListener interface
    virtual void handle(const Incident& inc);    

    /// Query interface - required of all Gaudi services
    virtual StatusCode queryInterface( const IID& riid, void** ppvUnknown );

    /// Provide the named ntuple Ptr from the data store -- dummy since we are not using this
    virtual SmartDataPtr<NTuple::Tuple> getNTuple(const char *tupleName){
        INTupleSvc *ntupleSvc=0;
        return SmartDataPtr<NTuple::Tuple> (ntupleSvc, "");
    }

    /// add a new item to an ntuple -- not supported
    virtual StatusCode addItem(const char *tupleName, 
        const char *item, double val){return StatusCode::FAILURE;};

    /** @brief Adds a <EM>pointer</EM> to an item -- only way to fill this guy
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column
    @param pval - pointer to a double value
    */
    virtual StatusCode addItem(const std::string & tupleName, 
        const std::string& itemName, const double* val);

    /// force writing of the ntuple to disk -- not suppored, but harmless since it happens anyway
    virtual StatusCode saveNTuples(){return StatusCode::SUCCESS;};

    /// Set a flag to denote whether or not to store a row at the end of this event
    virtual void storeRowFlag(bool flag) { m_storeFlag = flag; };
    /// retrieve the flag that denotes whether or not to store a row
    virtual bool storeRowFlag() { return m_storeFlag; };

private:
    /// Allow only SvcFactory to instantiate the service.
    friend class SvcFactory<RootTupleSvc>;

    RootTupleSvc ( const std::string& name, ISvcLocator* al );    

    bool m_storeFlag;

    /// routine to be called at the beginning of an event
    void beginEvent();
    /// routine that is called when we reach the end of an event
    void endEvent();

    StringProperty m_filename;
    StringProperty m_treename;
    StringProperty m_title;

    /// the ROOT stuff: a file and a a set of trees to put into it
    TFile * m_tf;

    std::map<std::string, TTree *> m_tree;

    int m_trials; /// total number of calls
    bool m_defaultStoreFlag;

    std::vector<float> m_floats; // needed for communication with ROOT's float branch.
    std::vector<const double*> m_pdoubles; // ordered  list of pointers to user variables

};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// declare the service factories for the ntupleWriterSvc
static SvcFactory<RootTupleSvc> a_factory;
const ISvcFactory& RootTupleSvcFactory = a_factory;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//         Implementation of RootTupleSvc methods
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RootTupleSvc::RootTupleSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc), m_trials(0)
{
    // declare the properties and set defaults
    declareProperty("filename",  m_filename="RootTupleSvc.root");
    declareProperty("treename", m_treename="1");
    declareProperty("title", m_title="Glast tuple");
    declareProperty("defaultStoreFlag", m_defaultStoreFlag=false);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::initialize () 
{
    StatusCode  status =  Service::initialize ();

    // bind all of the properties for this service
    setProperties ();
    std::string filename(m_filename);
    facilities::Util::expandEnvVar(&filename);

    // open the message log
    MsgStream log( msgSvc(), name() );

    // use the incident service to register "end" events
    IIncidentSvc* incsvc = 0;
    status = service("IncidentSvc", incsvc, true);

    if( status.isFailure() ) return status;

    incsvc->addListener(this, "BeginEvent", 100);
    incsvc->addListener(this, "EndEvent", 0);

    // -- set up the tuple ---
    m_tf   = new TFile( m_filename.value().c_str(), "RECREATE");
    // with the default treename, and default title
    m_tree[m_treename.value().c_str()] = new TTree( m_treename.value().c_str(),  m_title.value().c_str() );

    m_floats.reserve(250); // a little weakness: if larger than this, trouble
    return status;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::addItem(const std::string & tupleName, 
                                 const std::string& itemName, const double* val)
{
     MsgStream log(msgSvc(),name());
    StatusCode status = StatusCode::SUCCESS;
    m_floats.push_back(0);
    std::string treename=tupleName.empty()? m_treename.value() : tupleName;
    if( m_tree.find(treename)==m_tree.end()){
        // create new tree
        m_tree[treename]=new TTree(treename.c_str(), m_title.value().c_str());
        log << "Creating new tree " << treename << endreq;
    }
    m_tree[treename]->Branch(itemName.c_str(), &m_floats.back(), itemName.c_str());
    m_pdoubles.push_back(val);
    return status;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void RootTupleSvc::handle(const Incident &inc)
{
    // Purpose and Method:  This routine is called when an "incident"
    //   occurs.  This method determines what action the RootTupleSvc
    //   will take in response to a particular event.  Currently, we handle
    //   BeginEvent and EndEvent events.

    if(inc.type()=="BeginEvent") beginEvent();
    if(inc.type()=="EndEvent") endEvent();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void RootTupleSvc::beginEvent()
{
    /// Assume that we will NOT write out the row
    storeRowFlag(m_defaultStoreFlag);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void RootTupleSvc::endEvent()  // must be called at the end of an event to update, allow pause
{         
    ++m_trials;
    if (m_storeFlag == false) return;
    std::vector<float>::iterator fit = m_floats.begin();
    for( std::vector<const double*>::const_iterator pit = m_pdoubles.begin(); pit!=m_pdoubles.end(); ++pit){
        double t = static_cast<float>(**pit);
        *fit++ = static_cast<float>(**pit);
    }

    for( std::map<std::string, TTree*>::iterator it = m_tree.begin(); it!=m_tree.end(); ++it){
        TTree* t = it->second;
        t->Fill();
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::queryInterface(const IID& riid, void** ppvInterface)  {
    if ( IID_INTupleWriterSvc.versionMatch(riid) )  {
        *ppvInterface = (INTupleWriterSvc*)this;
    }
    else  {
        return Service::queryInterface(riid, ppvInterface);
    }
    addRef();
    return SUCCESS;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

StatusCode RootTupleSvc::finalize ()
{
    // open the message log
    MsgStream log( msgSvc(), name() );

    for( std::map<std::string, TTree*>::iterator it = m_tree.begin(); it!=m_tree.end(); ++it){
        TTree* t = it->second; 

        if( t->GetEntries() ==0 ) {

            log << MSG::INFO << "No entries added to the TTree \"" << it->first <<"\" : not writing it" << endreq;

        }else{
            log << MSG::INFO << "Writing the TTree \"" << it->first<< "\" in file "<<m_filename.value() 
                << " with " 
                << t->GetEntries() << " rows (" << m_trials << " total events)"<< endreq;
            t->Print(); // make a summary
            t->Write();
        }
    }

    TDirectory *saveDir = gDirectory; 
    m_tf->cd(); 
    m_tf->Write(0,TObject::kOverwrite); 
    m_tf->Close(); 
    saveDir->cd();


    return StatusCode::SUCCESS;;
}