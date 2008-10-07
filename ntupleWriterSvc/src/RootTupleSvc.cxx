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
#include "GaudiKernel/MsgStream.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"
#include "facilities/Util.h"

// root includes
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
#include "TLeafD.h"
#include "TLeaf.h"

#include <map>
#include <fstream>
#include <iomanip>
#include <list>
#include <string>
#include <utility>

#ifdef WIN32
#include <float.h> // used to check for NaN
#else
#include <cmath>
#endif

namespace {

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
    virtual StatusCode queryInterface( const InterfaceID& riid, void** ppvUnknown );

// LSR 14-Jul-08 code for ntuple types

    /** @brief Adds a pointer to a double, or an array of doubles
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column. append [n] to make a fixed array of length n
    @param pval - pointer to a double value
    @param fileName - name of ROOT file: if it does not exist, it will be created
    */
    virtual StatusCode addItem(const std::string & tupleName, 
        const std::string& itemName, const double* pval, 
        const std::string& fileName=std::string(""));

    /** @brief Adds a pointer to a float
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column. append [n] to make a fixed array of length n
    @param pval - pointer to a float value
    @param fileName - name of ROOT file: if it does not exist, it will be created
    */
    virtual StatusCode addItem(const std::string & tupleName, 
        const std::string& itemName, const float* pval,
        const std::string& fileName=std::string(""));

    /** @brief Adds a pointer to an int 
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column. append [n] to make a fixed array of length n
    @param pval - pointer to a int value
    @param fileName - name of ROOT file: if it does not exist, it will be created
    */
    virtual StatusCode addItem(const std::string & tupleName, 
        const std::string& itemName, const int* pval,
        const std::string& fileName=std::string(""));

    /** @brief Adds a pointer to an unsigned int 
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column. append [n] to make a fixed array of length n
    @param pval - pointer to a int value
    @param fileName - name of ROOT file: if it does not exist, it will be created
    */
    virtual StatusCode addItem(const std::string & tupleName, 
        const std::string& itemName, const unsigned int* pval,
        const std::string& fileName=std::string(""));

    /** @brief Adds a pointer to an unsigned long long 
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column. append [n] to make a fixed array of length n
    @param pval - pointer to a int value
    @param fileName - name of ROOT file: if it does not exist, it will be created
    */
    virtual StatusCode addItem(const std::string & tupleName, 
        const std::string& itemName, const unsigned long long* pval,
        const std::string& fileName=std::string(""));

    /** @brief Adds a pointer to a zero-terminated array of char 
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column.
    @param pval - pointer to the character array
    @param fileName - name of ROOT file: if it does not exist, it will be created
    */
    virtual StatusCode addItem(const std::string & tupleName, 
        const std::string& itemName, const char * val,
            const std::string& fileName=std::string(""));



    /** @brief interface to ROOT to add any item
    @param tupleName - name of the Root tree: if it does not exist, it will be created. If blank, use the default
    @param itemName - name of the tuple column. append [n] to make a fixed array of length n
    @param type  - ROOT-specific type qualifier ("/F", for example)
    @param pval - pointer to a void value or array
    @param fileName - name of ROOT file: if it does not exist, it will be created
    */
    StatusCode addAnyItem(const std::string & tupleName, 
               const std::string& itemName, const std::string type, 
               const void* pval, const std::string& fileName=std::string(""));

    /// Set a flag to denote whether or not to store a row at the end of this event,
    virtual void storeRowFlag(bool flag) { m_storeAll = flag; }
    /// retrieve the flag that denotes whether or not to store a row
    virtual bool storeRowFlag() { return m_storeAll; }

    // Had to remove const due to reprocessing needs, and the requirement to store branch pointers in some cases
    virtual std::string getItem(const std::string & tupleName, 
        const std::string& itemName, void*& pval);

    /** store row flag by tuple Name option, retrive currrent
    @param tupleName Name of the tuple (TTree for RootTupleSvc implemetation)
    @param flag new value
    @return previous value
    If service does not implement, it is ignored (return false)
    */
    virtual bool storeRowFlag(const std::string& tupleName, bool flag);


    //! Returns a pointer to the requested input TTree
    virtual bool getInputTreePtr(void*& treePtr,
                         const std::string& tupleName="MeritTuple");

    //! Returns a pointer to the requested output TTree
    virtual bool getOutputTreePtr(void*& treePtr,
                         const std::string& tupleName="MeritTuple");

    //! Save the row in the output file
    virtual void saveRow(const std::string& tupleName);

    /// allow clients to set TTree buffer size on a per branch basis or
    /// for whole TTree by setting bname="*"
    virtual void setBufferSize(const std::string& tupleName, int bufSize, 
                               const std::string& bname=std::string("*"));

private:
    /// Allow only SvcFactory to instantiate the service.
    friend class SvcFactory<RootTupleSvc>;

    RootTupleSvc ( const std::string& name, ISvcLocator* al );    

    StatusCode checkForNAN(TTree*, MsgStream& log);

    bool fileExists( const std::string & filename );

    /// For getting "the" current tree...
    bool getTree(std::string& treeName, TTree*& t);

    /// routine to be called at the beginning of an event
    void beginEvent();
    /// routine that is called when we reach the end of an event
    StatusCode endEvent();

    // Associated with the name of the first output ROOT file
    StringProperty m_filename;
    StringArrayProperty m_inFileJoParam;
    // stores the list of input files after env variables have been expanded
    std::vector<std::string> m_inFileList;
    StringProperty m_treename;
    StringProperty m_title;

    StringProperty m_jobInfoTreeName;
    StringProperty m_jobInfo;

    /// the ROOT stuff: a file and a a set of trees to put into it
    // replaced with m_fileCol, so we can handle multiple ROOT output files
    std::map<std::string, TFile*> m_fileCol;

    /// collection of output TTrees
    std::map<std::string, TTree *> m_tree;

    //std::map<std::string, TTree *> m_inTree;

    /// collection of input TChains
    std::map<std::string, TChain *> m_inChain;

    /// collection of leaf addresses for items that we have to create an object for.
    /// This occurs in reprocessing, when not all AnaTup Tools are executed - so not all branches have a corresponding
    /// variable to TChain::SetBranchAddress for, so that we have a stable location to provide via the getItem call.
    std::map<std::string, void*> m_itemPool;

    /// the flags, one per tree, for storing at the end of an event
    // assumes each TTree has a unique name
    std::map<std::string, bool> m_storeTree;

    /// If reading an input tuple also, then this is next event
    long long m_nextEvent;

    /// if set, store all ttrees 
    bool m_storeAll;

    int m_trials; /// total number of calls
    bool m_defaultStoreFlag;
    IntegerProperty m_autoSave; // passed to TTree::SetAutoSave.

    /// keep track of how many events had non-finite values
    int m_badEventCount;
    // assumes each leaf has a uniue name across all trees and files
    std::map<std::string, int> m_badMap; ///< map of counts for individual values
    BooleanProperty m_rejectIfBad; ///< set true to reject the tuple entry if bad values

    /// JO parameter to set the default buffer size for all TTrees
    int m_bufferSize;

    Char_t test; // to see what this is

    int m_updateCount;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// declare the service factories for the ntupleWriterSvc
static SvcFactory<RootTupleSvc> a_factory;
const ISvcFactory& RootTupleSvcFactory = a_factory;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//         Implementation of RootTupleSvc methods
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RootTupleSvc::RootTupleSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc), m_nextEvent(0), m_trials(0), 
  m_badEventCount(0)
{
    // declare the properties and set defaults
    declareProperty("filename",  m_filename="RootTupleSvc.root");

    // Provide a default empty list for inFileList parameter
    StringArrayProperty initList;
    std::vector<std::string> initVec;
    initVec.clear();
    initList.setValue(initVec);
    declareProperty("inFileList",  m_inFileJoParam=initList);

    declareProperty("treename", m_treename="1");
    declareProperty("title", m_title="Glast tuple");
    declareProperty("defaultStoreFlag", m_defaultStoreFlag=false);
    declareProperty("AutoSave", m_autoSave=100000); // ROOT default is 10000000
    declareProperty("RejectIfBad", m_rejectIfBad=true); 
    declareProperty("JobInfoTreeName", m_jobInfoTreeName="jobinfo");
    declareProperty("JobInfo", m_jobInfo=""); // string, if present, will write out single TTree entry
    declareProperty("BufferSize",m_bufferSize=32000);
    declareProperty("StartingIndex",m_nextEvent=0);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::initialize () 
{
    StatusCode  status =  Service::initialize ();

    // shut down ROOT's exception handling for Gleam runs
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

    m_fileCol.clear();
    m_tree.clear();
    //m_inTree.clear();
    m_badMap.clear();
    m_inChain.clear();
    m_inFileList.clear();
    m_itemPool.clear();

    // Split here depending on whether we are reading an input ntuple
    // and augmenting its output
    TDirectory* curdir = gDirectory; // will prevent unauthorized use

    /* HMK Not adding input TFiles to the m_fileCol, since they will 
       be apart of the TChain.
       We will expand any env variables in the input JO parameter.
       */
    if (m_inFileJoParam.value().size() > 0)
    {
        facilities::Util::expandEnvVarList(m_inFileJoParam.value(),m_inFileList);
        std::vector<std::string>::iterator it;
        for(it = m_inFileList.begin(); it != m_inFileList.end(); it++) {
            std::string tempStr = *it;
            facilities::Util::expandEnvVar(&tempStr);
            *it = tempStr;
        }
       /*
        TFile* tf = new TFile(m_inFileList.value().c_str(), "READ");
        if (!tf->IsOpen())
        {
            log << MSG::ERROR 
                << "cannot open ROOT file: " << m_filename.value() << endreq;
            delete tf;
            return StatusCode::FAILURE;
        }
        m_fileCol[m_inFilename.value()] = tf;
        */
    }

    // -- create primary output root file---
    TFile *tf   = new TFile( m_filename.value().c_str(), "RECREATE");
    if (!tf->IsOpen()) {
        log << MSG::ERROR 
            << "cannot open ROOT file: " << m_filename.value() << endreq;
        delete tf;
        return StatusCode::FAILURE;
    }
    m_fileCol[m_filename.value()] = tf;

    curdir->cd(); // restore previous directory

    m_updateCount = 0;

    return status;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bool RootTupleSvc::getTree(std::string& treeName, TTree*& t)
{
    MsgStream log(msgSvc(),name());
    TDirectory* saveDir = gDirectory; // will prevent unauthorized use

    bool inFileFlag = false;
    t = 0;

    // Are we reading from input ntuple(s)?
    if (m_inFileList.size() > 0)
    {
        inFileFlag = true;
        std::map<std::string, TChain*>::iterator inIter = m_inChain.find(treeName);
        // Must self-initialize the first time around
        if (inIter == m_inChain.end()) {
            
            // Assuming all input ROOT files have the same tree name,
            // we could adjust the code to pick out the TTree and provide 
            // its name in the TChain::Add call
            TChain *ch = new TChain(treeName.c_str());
            // Add all files in the JO list
            std::vector<std::string>::const_iterator fileListItr;
            for (fileListItr = m_inFileList.begin();
                 fileListItr != m_inFileList.end();
                 fileListItr++ ) {
                std::string fileName = *fileListItr;
                facilities::Util::expandEnvVar(&fileName);
                    int stat = ch->Add(fileName.c_str());
                    if (stat <= 0)
                        log << MSG::WARNING << "Failed to TChain::Add " 
                            << fileName << " return code: " << stat << endreq;
                    else
                        log << MSG::INFO << "Added File: " << fileName << endreq;

            }  // end for loop TChain initialized
            // add new TChain to the map
            m_inChain[treeName] = ch;
            // call GetEntries to load the headers of the TFiles
            long long nevents = ch->GetEntries();
            log << MSG::INFO << "Number of events in input files = " 
                << nevents << " StartingIndex: " << m_nextEvent << endreq;
            if ((m_nextEvent > nevents-1) || (m_nextEvent < 0)) {
                log << MSG::WARNING << "StartingIndex invalid, resetting "
                    << m_nextEvent << " to zero" << endreq;
                m_nextEvent = 0;
            }
            int numbytes = ch->GetEntry(m_nextEvent);
            if (numbytes <= 0) 
                log << MSG::WARNING << "Unable to read tuple event, "
                    << m_nextEvent << endreq;
            // Here is our chance to set up branch pointers for the whole input TChain, so that
            // no elements are missed
            TObjArray *brCol = ch->GetListOfBranches();
            int numBranches = brCol->GetEntries();
            int iBranch;
            for (iBranch=0;iBranch<numBranches;iBranch++) {
                std::string branchName(((TBranch*)(brCol->At(iBranch)))->GetName());
                std::string leafName = branchName;
                int index = leafName.find("[");
                if(index!=std::string::npos) leafName = leafName.substr(0,index); 
                log << MSG::DEBUG << "setting branch: " << branchName
                    << " and Leaf: " << leafName << endreq;
                TLeaf *leaf = ((TBranch*)brCol->At(iBranch))->GetLeaf(leafName.c_str());
                if (!leaf) {
                    log << MSG::WARNING << "Leaf: " << leafName << " not found" << endreq;
                    continue;
                }
                std::string type_name = leaf->GetTypeName();
                if (type_name == "Float_t")
                {
                    m_itemPool[branchName] = new Float_t[leaf->GetNdata()];
                }
                else if (type_name == "Int_t")
                {
                    m_itemPool[branchName]= new Int_t[leaf->GetNdata()];
                }
                else if (type_name == "UInt_t")
                {
                    m_itemPool[branchName]= new UInt_t[leaf->GetNdata()];
                }        
                else if (type_name == "ULong64_t")
                {
                    m_itemPool[branchName]= new ULong64_t[leaf->GetNdata()];
                }
                else if (type_name == "Double_t")
                {
                    m_itemPool[branchName] = new Double_t[leaf->GetNdata()];
                }
                else if (type_name == "Char_t")
                {
                    m_itemPool[branchName] = new Char_t[leaf->GetNdata()];
                }
                else
                {
                    log << MSG::WARNING << "type: " << type_name <<" not found" << endreq;
                }
                ch->SetBranchAddress(branchName.c_str(), m_itemPool[branchName]);
            }

            inIter = m_inChain.find(treeName);

        } // end if for initialization first time 
        // copy the current tree to allow access to its branches
        // zero is the number of entries to copy - so we're just 
        // copying the TTree structure not contents
        t=inIter->second->CloneTree(0);
    } // end check for input files

    // Back to regular situation, where we are setting up an output file
    if (t == 0)
    {
        t = new TTree(treeName.c_str(), m_title.value().c_str());
        t->SetAutoSave(m_autoSave);
    }

    saveDir->cd();
    return inFileFlag;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

StatusCode RootTupleSvc::addAnyItem(const std::string & tupleName, 
                                    const std::string& itemName0, 
                                    const std::string type,  
                                    const void* pval, 
                                    const std::string& fileName)
{
    MsgStream log(msgSvc(),name());
    //bool inFileFlag = false;
    StatusCode status = StatusCode::SUCCESS;
    std::string treename=tupleName.empty()? m_treename.value() : tupleName;
    std::string rootFileName = fileName.empty() ? m_filename.value() : fileName;
    TDirectory *saveDir = gDirectory;

    // Check list of output files
    if ( m_fileCol.find(rootFileName) == m_fileCol.end()) {
        // create a new TFile
        TFile *tf = new TFile(rootFileName.c_str(), "RECREATE");
        if (!tf->IsOpen()) {
            log << MSG::ERROR 
                << "cannot open ROOT file: " << rootFileName << endreq;
            delete tf;
            saveDir->cd();
            return StatusCode::FAILURE;
        }
        m_fileCol[rootFileName] = tf;
        getTree(treename,m_tree[treename]);
        m_tree[treename]->SetDirectory(tf);
        log << MSG::INFO << "Creating new tree \"" << treename << "\"" 
            << " in file: " << rootFileName << endreq;
    } else if( m_tree.find(treename)==m_tree.end()){
        // create new tree
        m_fileCol[rootFileName]->cd();
        getTree(treename,m_tree[treename]);
        m_tree[treename]->SetDirectory(m_fileCol[rootFileName]);
        log << MSG::INFO << "Creating new tree \"" << treename << "\"" 
            << " in file: " << rootFileName << endreq;
    }
    // note have to cast away the const here!

    // Searches list of branches, and returns NULL if itemName0 is not found
    TBranch* thisBranch = m_tree[treename]->GetBranch(itemName0.c_str());
    if(thisBranch==NULL) {
        log << MSG::DEBUG << "Creating new branch in AddAny for " << itemName0
            << endreq;
        // This is a new branch
        m_tree[treename]->Branch(itemName0.c_str(), 
            const_cast<void*>(pval), (itemName0+type).c_str(),m_bufferSize);
    } else {
        log << MSG::DEBUG << "Found branch in TTree: " << itemName0
            << endreq;
        thisBranch->SetAddress(const_cast<void*>(pval));
    }
    saveDir->cd();
    return status;
}
/* about these codes:

            - C : a character string terminated by the 0 character
            - B : an 8 bit signed integer (Char_t)
            - b : an 8 bit unsigned integer (UChar_t)
            - S : a 16 bit signed integer (Short_t)
            - s : a 16 bit unsigned integer (UShort_t)
            - I : a 32 bit signed integer (Int_t)
            - i : a 32 bit unsigned integer (UInt_t)
            - F : a 32 bit floating point (Float_t)
            - D : a 64 bit floating point (Double_t)
            - L : a 64 bit signed integer (Long64_t)
            - l : a 64 bit unsigned integer (ULong64_t)
            - O : a boolean (Bool_t)

  */

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::addItem(const std::string & tupleName, 
                                 const std::string& itemName, 
                                 const double* pval, 
                                 const std::string& fileName)
{
    return addAnyItem(tupleName, itemName, "/D", (void*)pval, fileName);


}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::addItem(const std::string & tupleName, 
                                 const std::string& itemName, 
                                 const float* pval, 
                                 const std::string& fileName)
{
    return addAnyItem(tupleName, itemName, "/F", (void*)pval, fileName);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::addItem(const std::string & tupleName, 
                                 const std::string& itemName, 
                                 const int* pval, 
                                 const std::string& fileName)
{
    return addAnyItem(tupleName, itemName, "/I", (void*)pval, fileName);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::addItem(const std::string & tupleName, 
                                 const std::string& itemName, 
                                 const unsigned int* pval, 
                                 const std::string& fileName)
{
    return addAnyItem(tupleName, itemName, "/i", (void*)pval, fileName);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::addItem(const std::string & tupleName, 
                                 const std::string& itemName, 
                                 const unsigned long long* pval, 
                                 const std::string& fileName)
{
    return addAnyItem(tupleName, itemName, "/l", (void*)pval, fileName);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::addItem(const std::string & tupleName, 
                               const std::string& itemName,
                               const char * pval,
                               const std::string& fileName)
{
    return addAnyItem(tupleName, itemName, "/C", (void*)pval, fileName);
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
    TDirectory *saveDir = gDirectory;
    /// If we have an input ntuple then read the branches...
    for(std::map<std::string, TChain*>::iterator inIter = m_inChain.begin(); inIter != m_inChain.end(); inIter++)
    {
        MsgStream log(msgSvc(),name());
        std::string treeName = inIter->first;
        //inIter->second->LoadTree(m_nextEvent);
        log << MSG::DEBUG << "event: " << m_nextEvent << " TreeNum: "
            << inIter->second->GetTreeNumber() << endreq;

        
        int numBytes = inIter->second->GetEntry(m_nextEvent++);
        if (numBytes <= 0){
            MsgStream log(msgSvc(),name());
            log << MSG::ERROR << "Failed to load event " << m_nextEvent-1
                << " from the input chain, terminating job" << endreq;
            exit(1);
        }

    }

    /// Assume that we will NOT write out the row
    storeRowFlag(m_defaultStoreFlag);
    for(std::map<std::string, bool>::iterator it=m_storeTree.begin(); it!=m_storeTree.end(); ++it){
        it->second=false;
    }

    saveDir->cd();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::endEvent()
    // must be called at the end of an event to update, allow pause
{         
    MsgStream log(msgSvc(),name());
    StatusCode sc = SUCCESS;
    TDirectory *saveDir = gDirectory;

    ++m_trials;
    for( std::map<std::string, TTree*>::iterator it = m_tree.begin();
         it!=m_tree.end(); ++it){
        if( m_storeAll || m_storeTree[it->first]  ) {
            TTree* t = it->second;
            t->GetCurrentFile()->cd();
            sc = checkForNAN(t, log);
            // check the tuple for non-finite entries, do not fill the tuple if found (unless overriden)
            if( sc.isFailure() ){ 
                m_badEventCount++; 
                if (!m_rejectIfBad) t->Fill();
            }else{
                t->Fill();
            }
            m_storeTree[it->first]=false;
        }
    }
        
    saveDir->cd();
    return sc;

}
StatusCode RootTupleSvc::checkForNAN( TTree* t, MsgStream& log)
{
    TDirectory *saveDir = gDirectory;
    t->GetCurrentFile()->cd();
    TObjArray* ta = t->GetListOfBranches();
    StatusCode sc = SUCCESS;

    // now iterate.
    int entries = ta->GetEntries();
    for( int i = 0; i<entries; ++i) { // should try a TIter
        TBranch * b = (TBranch*)(*ta)[i];
        TLeaf* leaf = (TLeaf*)(*b->GetListOfLeaves())[0]; 
        double val = leaf->GetValue();
        log << MSG::DEBUG << leaf->GetName() << " val: " << val << endreq;
        void* valpointer = leaf->GetValuePointer();
        if( ! isFinite(val) ){
            log << MSG::DEBUG  << "Tuple item " << leaf->GetName() << " is not finite!" << endreq;
            m_badMap[leaf->GetName()]++;
            sc = StatusCode::FAILURE;
        }
    }
    saveDir->cd();
    return sc;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode RootTupleSvc::queryInterface(const InterfaceID& riid, void** ppvInterface)  {
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

    // -- set up job info TTree if requested to add values, or the tree exists already

    TTree * jobinfotree(0);
    std::map<std::string, TTree*>::const_iterator treeit = m_tree.find(m_jobInfoTreeName);
    if( treeit!=m_tree.end()){
        jobinfotree= treeit->second;
    }

    std::list<float> values;
    if( jobinfotree!=0 || ! m_jobInfo.value().empty() ){
        std::map<std::string, std::string > parmap;
        facilities::Util::keyValueTokenize(m_jobInfo.value(), ",", parmap);
        for(  std::map<std::string, std::string >::const_iterator mip = parmap.begin(); mip!=parmap.end(); ++mip){
            try {
            std::string key (mip->first), value(mip->second);
            // assume all numbers
            values.push_back( facilities::Util::stringToDouble(value) );
            float* fval = &values.back(); 
            addItem(m_jobInfoTreeName, key, fval);
            } catch(...) {
               log << MSG::WARNING << "RootTupleSvc exception caught "
                   << "processing jobinfo " << m_jobInfo.value()  
                   << "Check string is in form: a=1,b=2,c=3" <<endreq; 
            }
        }
        // process this row
        saveRow(m_jobInfoTreeName); 
        log << MSG::INFO << "jobinfo scan: "<< endreq;
        if( jobinfotree!=0) jobinfotree->Scan();
    }

    for( std::map<std::string, TTree*>::iterator it = m_tree.begin(); it!=m_tree.end(); ++it){
        TTree* t = it->second; 
        t->GetCurrentFile()->cd();
        if( m_storeTree[it->first] ) t->Fill(); // In case the algorithm did an entry during its finalize

        if( t->GetEntries() ==0 ) {

            log << MSG::INFO << "No entries added to the TTree \"" << it->first <<"\" : not writing it" << endreq;

        }else{
            log << MSG::INFO << "Writing the TTree \"" << it->first<< "\" in file "<< t->GetCurrentFile()->GetName() //m_filename.value() 
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
        log << MSG::WARNING << "==================================================================" << endreq;
        log << MSG::WARNING << "Found " << m_badEventCount << " bad events: table of bad values follows\n" ;
        if( log.isActive()) {
            log.stream() << "\t\t\t" << std::setw(20)<< std::left << "Name" 
                << std::setw(10) << std::right << "count" << std::endl;
            for( std::map<std::string,int>::const_iterator it=m_badMap.begin(); it !=m_badMap.end(); ++it){
                log.stream() << "\t\t\t"<< std::setw(20) << std::left << it->first 
                    << std::setw(10)<<  std::right << it->second << std::endl;
            }
        }
        log << endreq;
        if (m_rejectIfBad) { 
            log << MSG::WARNING << "==========> REJECTED since RejectIfBad flag set to do so!" << endreq;
        } 
    }

 


    for( std::map<std::string, TFile*>::iterator it = m_fileCol.begin(); it!=m_fileCol.end(); ++it){
        TFile* f = it->second; 
        if ((!f) || (!f->IsOpen())) {
            log << MSG::WARNING << "ROOT File: " << f->GetName() 
                << " is not open - skipping write" << endreq;
        } else {
            f->cd();
            f->Write(0,TObject::kOverwrite);
            f->Close();
        }
        if (f) {
            delete f;
            f=0;
        }
    }

    return StatusCode::SUCCESS;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool RootTupleSvc::storeRowFlag(const std::string& tupleName, bool flag)
{
    bool t = m_storeTree[tupleName];
    m_storeTree[tupleName] = flag;
    return t;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool RootTupleSvc::getInputTreePtr(void*& pval, const std::string & tupleName)
{
    MsgStream log(msgSvc(),name());

    std::string treename=tupleName.empty()? m_treename.value() : tupleName;
    TDirectory *saveDir = gDirectory;

    // Check input files for this TTree
    std::map<std::string, TChain*>::const_iterator chainit = m_inChain.find(treename);
    if (chainit != m_inChain.end()) {
        // Found the TChain, now return
        TChain* ch = chainit->second;
        ch->GetCurrentFile()->cd();

        pval = (void *)(ch->GetTree());
        saveDir->cd();
        return true;

    }

    log << MSG::INFO << "Did not find tree" << treename << endreq;
    pval = 0;
    saveDir->cd();
    return false;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool RootTupleSvc::getOutputTreePtr(void*& pval, const std::string & tupleName)
{
    MsgStream log(msgSvc(),name());

    std::string treename=tupleName.empty()? m_treename.value() : tupleName;
    TDirectory *saveDir = gDirectory;

    // Check input files for this TTree
    std::map<std::string, TTree*>::const_iterator treeit = m_tree.find(treename);
    if (treeit != m_tree.end()) {
        // Found the TChain, now return
        TTree* t = treeit->second;
        t->GetCurrentFile()->cd();

        pval = (void *)t;
        saveDir->cd();
        return true;

    }

    log << MSG::INFO << "Did not find tree" << treename << endreq;
    pval = 0;
    saveDir->cd();
    return false;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
std::string RootTupleSvc::getItem(const std::string & tupleName, 
                                   const std::string& itemName, void*& pval)
{
    // Hack for inputChain required for reprocessing option

    MsgStream log(msgSvc(),name());

    std::string treename=tupleName.empty()? m_treename.value() : tupleName;
    TDirectory *saveDir = gDirectory;

    std::map<std::string, TTree*>::const_iterator treeit = m_tree.find(treename);
    if( treeit==m_tree.end()){
        log << MSG::ERROR << "Did not find tree" << treename << endreq;
        throw std::invalid_argument("RootTupleSvc::getItem: did not find tuple or leaf");
    }
    TTree* t = treeit->second;
    t->GetCurrentFile()->cd();

    if( itemName.empty()){
        // assume this is a request for the tree
        pval = (void *)t;
        saveDir->cd();
        return itemName;
    }

    TLeaf *leaf = 0;
    bool foundInChain = false;
    // Check potential input tree 
    std::map<std::string, TChain*>::const_iterator inputChain = m_inChain.find(treename);
    if (inputChain != m_inChain.end()) 
        leaf = inputChain->second->GetLeaf(itemName.c_str());
    

    if (leaf == 0) 
        leaf = t->GetLeaf(itemName.c_str());
    else
        foundInChain = true;

    if( leaf==0)
        throw std::invalid_argument(std::string("RootTupleSvc::getItem: did not find tuple or leaf: ")+ itemName);
    
    pval = leaf->GetValuePointer();
 
    std::string type_name(leaf->GetTypeName());

    if (foundInChain) {
        std::map<std::string, void*>::iterator itemIt = m_itemPool.find(itemName);
        // Create a new object to store this leaf pointer
        // This is necessary when we move to a new TTree in the TChain, otherwise, this address will be lost
        // and unusable by the clients that are relying on a stable address
        log << MSG::DEBUG << "item: " << itemName << " type: " << type_name << " dim: " << leaf->GetNdata() << endreq;
        if (itemIt == m_itemPool.end()) {
            if (type_name == "Float_t")
            {
                m_itemPool[itemName] = new Float_t[leaf->GetNdata()];
            }
            else if (type_name == "Int_t")
            {
                m_itemPool[itemName]= new Int_t[leaf->GetNdata()];
            }
            else if (type_name == "UInt_t")
            {
                m_itemPool[itemName]= new UInt_t[leaf->GetNdata()];
            }        
            else if (type_name == "ULong64_t")
            {
                m_itemPool[itemName]= new ULong64_t[leaf->GetNdata()];
            }
            else if (type_name == "Double_t")
            {
                m_itemPool[itemName] = new Double_t[leaf->GetNdata()];
            }
            else if (type_name == "Char_t")
            {
                m_itemPool[itemName] = new Char_t[leaf->GetNdata()];
            }
            else
            {
                log << MSG::WARNING << "type: " << type_name <<" not found" << endreq;
            }
            inputChain->second->SetBranchAddress(itemName.c_str(), m_itemPool[itemName]);
            leaf = inputChain->second->GetLeaf(itemName.c_str());
            pval = leaf->GetValuePointer();
        }
    }
    saveDir->cd();
    return type_name;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void RootTupleSvc::saveRow(const std::string& tupleName)
{
    std::map<std::string, TTree*>::iterator treeit=m_tree.find(tupleName);
    if( treeit==m_tree.end()){
        MsgStream log(msgSvc(),name());
        log << MSG::ERROR << "Did not find tree " << tupleName << endreq;
        throw std::invalid_argument("RootTupleSvc::saveRow: did not find tupleName");
    }

    TTree* t= treeit->second;
    t->Fill();
    m_storeTree[treeit->first]=false;
}

void RootTupleSvc::setBufferSize(const std::string& tupleName, int bufSize,
                                 const std::string& bname) {

    std::map<std::string, TTree*>::iterator treeit=m_tree.find(tupleName);
    if( treeit==m_tree.end()){
        MsgStream log(msgSvc(),name());
        log << MSG::ERROR << "Did not find tree " << tupleName << endreq;
        throw std::invalid_argument("RootTupleSvc::setBufferSize: did not find tupleName");
    }

    TTree* t= treeit->second;
    if (t)
        t->SetBasketSize(bname.c_str(), bufSize);
    else {
        MsgStream log(msgSvc(),name());
        log << MSG::WARNING << "Could not obtain tree, " << tupleName
                            << " Not setting buffer size" << endreq;
     }

}

bool RootTupleSvc::fileExists( const std::string & filename )
 {
  TDirectory *saveDir = gDirectory;
  bool fileExists = false ;
  TFile * file = TFile::Open(filename.c_str()) ;
  if (file)
   {
    if (!file->IsZombie()) 
     {
      file->Close() ;
      fileExists = true ;
     }
    delete file ;
   }
  saveDir->cd();
  return fileExists ;
 }  
