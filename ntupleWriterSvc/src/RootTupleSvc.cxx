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
#include "TSystem.h"


/** 
* @class RootTupleSvc
* @brief Special service that directly writes ROOT tuples
*
* It also allows multiple TTree's in the root file: see the addItem (by pointer) member function.
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
        const std::string& itemName, const double* pval);

    /// force writing of the ntuple to disk -- not suppored, but harmless since it happens anyway
    virtual StatusCode saveNTuples(){return StatusCode::SUCCESS;};

    /// Set a flag to denote whether or not to store a row at the end of this event,
    virtual void storeRowFlag(bool flag) { m_storeAll = flag; };
    /// retrieve the flag that denotes whether or not to store a row
    virtual bool storeRowFlag() { return m_storeAll; };

    /** store row flag by tuple Name option, retrive currrent
    @param tupleName Name of the tuple (TTree for RootTupleSvc implemetation)
    @param flag new value
    @return previous value
    If service does not implement, it is ignored (return false)
    */
    virtual bool storeRowFlag(const std::string& tupleName, bool flag);;

private:
    /// Allow only SvcFactory to instantiate the service.
    friend class SvcFactory<RootTupleSvc>;

    RootTupleSvc ( const std::string& name, ISvcLocator* al );    


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

    /// the flags, one per tree, for storing at the end of an event
    std::map<std::string, bool> m_storeTree;

    /// if set, store all ttrees 
    bool m_storeAll;

    int m_trials; /// total number of calls
    bool m_defaultStoreFlag;
    IntegerProperty m_autoSave; // passed to TTree::SetAutoSave.

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
    declareProperty("AutoSave", m_autoSave=100000); // ROOT default is 10000000


}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::initialize () 
{
    StatusCode  status =  Service::initialize ();

    gSystem->ResetSignal(kSigBus);
    gSystem->ResetSignal(kSigSegmentationViolation);
    gSystem->ResetSignal(kSigIllegalInstruction);
    gSystem->ResetSignal(kSigFloatingException); 

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
    //TTree* t = new TTree( m_treename.value().c_str(),  m_title.value().c_str() );
    //m_tree[m_treename.value().c_str()] = t;
    //t->SetAutoSave(m_autoSave); 

    return status;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::addItem(const std::string & tupleName, 
                                 const std::string& itemName, const double* pval)
{
    MsgStream log(msgSvc(),name());
    StatusCode status = StatusCode::SUCCESS;
    std::string treename=tupleName.empty()? m_treename.value() : tupleName;
    TDirectory *saveDir = gDirectory;
    if( m_tree.find(treename)==m_tree.end()){
        // create new tree
        m_tf->cd();
        TTree* t = new TTree(treename.c_str(), m_title.value().c_str());
        t->SetAutoSave(m_autoSave);
        m_tree[treename]=t;
        log << MSG::INFO << "Creating new tree \"" << treename << "\"" << endreq;
    }

    // this adds a branch with a pointer to a double (the "/D" after the second name)
    m_tree[treename]->Branch(itemName.c_str(), (void*)pval, (itemName+"/D").c_str());
    saveDir->cd();
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
    for(std::map<std::string, bool>::iterator it=m_storeTree.begin(); it!=m_storeTree.end(); ++it){
        it->second=false;
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void RootTupleSvc::endEvent()  // must be called at the end of an event to update, allow pause
{         
    ++m_trials;
    for( std::map<std::string, TTree*>::iterator it = m_tree.begin(); it!=m_tree.end(); ++it){
        if( m_storeAll || m_storeTree[it->first]  ) {
            TTree* t = it->second;
            t->Fill();
        }
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
            log << MSG::DEBUG;
            if( log.isActive() ){
                t->Print(); // make a summary (too bad ROOT doesn't allow you to specify a stream
            }
            log << endreq;
        }
    }

    TDirectory *saveDir = gDirectory; 
    m_tf->cd(); 
    m_tf->Write(0,TObject::kOverwrite); 
    m_tf->Close(); 
    saveDir->cd();


    return StatusCode::SUCCESS;;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bool RootTupleSvc::storeRowFlag(const std::string& tupleName, bool flag)
{
    bool t = m_storeTree[tupleName];
    m_storeTree[tupleName] = flag;
    return t;
}
