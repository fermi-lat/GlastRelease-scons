/** 
 * @file RootTupleSvc.cxx
 * @brief declare, implement the class RootTupleSvc
 *
 * Special service that directly writes ROOT tuples
 * It also allows multiple TTree's in the root file: see the addItem (by pointer) member function.
 * $Header$
 */

#include "GaudiKernel/Service.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"
#include "checkSum.h"
#include "facilities/Util.h"
#include <map>

// root includes
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TLeafD.h"
#include "TLeaf.h"

#include <fstream>

namespace {
#ifdef WIN32
#include <float.h> // used to check for NaN
#else
#include <cmath>
#endif

    bool isFinite(double val) {
        using namespace std; // should allow either std::isfinite or ::isfinite
#ifdef WIN32 
        return (_finite(val)!=0);  // Win32 call available in float.h
#else
        return (isfinite(val)!=0); // gcc call available in math.h 
#endif
    }
} // anom namespace




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

    /** @brief Adds a pointer to a double, or an array of doubles
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column. append [n] to make a fixed array of length n
    @param pval - pointer to a double value
    */
    virtual StatusCode addItem(const std::string & tupleName, 
        const std::string& itemName, const double* pval);

    /** @brief Adds a pointer to a float
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column. append [n] to make a fixed array of length n
    @param pval - pointer to a float value
    */
    virtual StatusCode addItem(const std::string & tupleName, 
        const std::string& itemName, const float* pval);

    /** @brief Adds a pointer to an int 
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column. append [n] to make a fixed array of length n
    @param pval - pointer to a int value
    */
    virtual StatusCode addItem(const std::string & tupleName, 
        const std::string& itemName, const int* pval);

    /** @brief interface to ROOT to add any item
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column. append [n] to make a fixed array of length n
    @param type  - ROOT-specific type qualifier ("/F", for example)
    @param pval - pointer to a void value or array
    */
    StatusCode addAnyItem(const std::string & tupleName, 
               const std::string& itemName, const std::string type, const void* pval);

    /// Set a flag to denote whether or not to store a row at the end of this event,
    virtual void storeRowFlag(bool flag) { m_storeAll = flag; }
    /// retrieve the flag that denotes whether or not to store a row
    virtual bool storeRowFlag() { return m_storeAll; }

    /** store row flag by tuple Name option, retrive currrent
    @param tupleName Name of the tuple (TTree for RootTupleSvc implemetation)
    @param flag new value
    @return previous value
    If service does not implement, it is ignored (return false)
    */
    virtual bool storeRowFlag(const std::string& tupleName, bool flag);

private:
    /// Allow only SvcFactory to instantiate the service.
    friend class SvcFactory<RootTupleSvc>;

    RootTupleSvc ( const std::string& name, ISvcLocator* al );    

    StatusCode checkForNAN(TTree*, MsgStream& log);

    /// routine to be called at the beginning of an event
    void beginEvent();
    /// routine that is called when we reach the end of an event
    StatusCode endEvent();

    StringProperty m_filename;
    StringProperty m_checksumfilename;
    StringProperty m_treename;
    StringProperty m_title;

    /// the ROOT stuff: a file and a a set of trees to put into it
    TFile * m_tf;

    /// the checksum object
    checkSum* m_checkSum;

    std::map<std::string, TTree *> m_tree;

    /// the flags, one per tree, for storing at the end of an event
    std::map<std::string, bool> m_storeTree;

    /// if set, store all ttrees 
    bool m_storeAll;

    int m_trials; /// total number of calls
    bool m_defaultStoreFlag;
    IntegerProperty m_autoSave; // passed to TTree::SetAutoSave.

    /// keep track of how many events rejected
    int m_badEventCount;

};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// declare the service factories for the ntupleWriterSvc
static SvcFactory<RootTupleSvc> a_factory;
const ISvcFactory& RootTupleSvcFactory = a_factory;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//         Implementation of RootTupleSvc methods
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RootTupleSvc::RootTupleSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc), m_trials(0), m_badEventCount(0)
{
    // declare the properties and set defaults
    declareProperty("filename",  m_filename="RootTupleSvc.root");
    declareProperty("checksumfilename", m_checksumfilename=""); // default empty
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
    if (!m_tf->IsOpen()) {
        log << MSG::ERROR 
            << "cannot open ROOT file: " << m_filename.value() << endreq;
        delete m_tf;
        return StatusCode::FAILURE;
    }
    // with the default treename, and default title
    //TTree* t = new TTree( m_treename.value().c_str(),  m_title.value().c_str() );
    //m_tree[m_treename.value().c_str()] = t;
    //t->SetAutoSave(m_autoSave); 

    // set up the check sum ofstream
    m_checkSum = new checkSum(m_checksumfilename);
    if ( m_checkSum->bad() ) {
        log << MSG::ERROR
            << "cannot open checksum file " << m_checksumfilename << endreq;
        delete m_checkSum;
        return StatusCode::FAILURE;
    }

    return status;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::addAnyItem(const std::string & tupleName, 
                                    const std::string& itemName, const std::string type,  const void* pval)
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
    // note have to cast away the const here!
    m_tree[treename]->Branch(itemName.c_str(), const_cast<void*>(pval), (itemName+type).c_str());
    saveDir->cd();
    return status;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::addItem(const std::string & tupleName, 
                                 const std::string& itemName, const double* pval)
{
    return addAnyItem(tupleName, itemName, "/D", (void*)pval);


}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::addItem(const std::string & tupleName, 
                                 const std::string& itemName, const float* pval)
{
    return addAnyItem(tupleName, itemName, "/F", (void*)pval);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::addItem(const std::string & tupleName, 
                                 const std::string& itemName, const int* pval)
{
    return addAnyItem(tupleName, itemName, "/I", (void*)pval);
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
StatusCode RootTupleSvc::endEvent()
    // must be called at the end of an event to update, allow pause
{         
    MsgStream log(msgSvc(),name());
    StatusCode sc = SUCCESS;

    ++m_trials;
    for( std::map<std::string, TTree*>::iterator it = m_tree.begin();
         it!=m_tree.end(); ++it){
        if( m_storeAll || m_storeTree[it->first]  ) {
            TTree* t = it->second;
            sc = checkForNAN(t, log);
            // check the tuple for non-finite entries, do not fill the tuple if found
            if( sc.isFailure() ){ 
                m_badEventCount++; 
            }else{
                t->Fill();
            }
            m_storeTree[it->first]=false;
            // doing the checksum here
            std::string treeName = t->GetName();
            if ( m_checkSum->is_open() && treeName == "MeritTuple" ) {
                log << MSG::VERBOSE << "calculating checksum for "
                    << treeName << endreq;
                m_checkSum->write(t);
            }
        }
    }
        
    return sc;

}
StatusCode RootTupleSvc::checkForNAN( TTree* t, MsgStream& log)
{
    TObjArray* ta = t->GetListOfBranches();
    StatusCode sc = SUCCESS;

    // now iterate.
    int entries = ta->GetEntries();
    for( int i = 0; i<entries; ++i) { // should try a TIter
        TBranch * b = (TBranch*)(*ta)[i];
        TLeaf* leaf = (TLeaf*)(*b->GetListOfLeaves())[0]; 
        double val = leaf->GetValue();
        if( ! isFinite(val) ){
            log << MSG::ERROR  << "Tuple item " << leaf->GetName() << " is not finite!" << endreq;
            sc = StatusCode::FAILURE;
        }
    }
    return sc;
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
        if( m_storeTree[it->first] ) t->Fill(); // In case the algorithm did an entry during its finalize

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
    if (m_badEventCount>0){
        log << MSG::ERROR << "Found and rejected " << m_badEventCount << " bad events!" << endreq;
    }

    TDirectory *saveDir = gDirectory; 
    m_tf->cd(); 
    m_tf->Write(0,TObject::kOverwrite); 
    m_tf->Close(); 
    saveDir->cd();

    // deleting the checksum object
    if ( m_checkSum )
        delete m_checkSum;

    return StatusCode::SUCCESS;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bool RootTupleSvc::storeRowFlag(const std::string& tupleName, bool flag)
{
    bool t = m_storeTree[tupleName];
    m_storeTree[tupleName] = flag;
    return t;
}
