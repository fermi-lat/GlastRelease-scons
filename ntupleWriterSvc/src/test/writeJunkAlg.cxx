#define ntupleWriter_writeJunkAlg_CPP 


#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "TTree.h"


#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/StatusCode.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"
#include <cmath>

/**
 * @class writeJunkAlg
 * @brief test algorithm for the ntupleWriterSvc 
 *
 * This algorithm tests the creation and writing of ntuples via
 * the ntupleWriterSvc.  The output from this routine is a
 * ROOT ntuple, containing 19 entries.
 */

class writeJunkAlg : public Algorithm {
    
public:
    /// Constructor of this form must be provided for all Gaudi algorithms
    writeJunkAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    std::string m_tupleName;
    float m_float;
    int   m_int;
    unsigned int m_uint;
    double m_count; // special value to put into tuples
    double m_square; // another
   
    float m_float2;

    double m_array[2]; // test an array

    float m_memoryFloat;
    int m_memoryInt;

    INTupleWriterSvc *m_rootTupleSvc;

    float* m_float_test;

    char  m_name[10];

};

//static const AlgFactory<writeJunkAlg>  Factory;
//const IAlgFactory& writeJunkAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(writeJunkAlg);

writeJunkAlg::writeJunkAlg(const std::string& name, ISvcLocator* pSvcLocator) 
: Algorithm(name, pSvcLocator)
{

    strncpy(m_name, "default",10); // make it point to something

    declareProperty("tupleName",  m_tupleName="");

}

StatusCode writeJunkAlg::initialize() {
    
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

      // get a pointer to RootTupleSvc as well
    sc = service("RootTupleSvc", m_rootTupleSvc);

    if( sc.isFailure() ) {
        log << MSG::ERROR << "writeJunkAlg failed to get the RootTupleSvc" << endreq;
        return sc;
    }

    m_count = 0;

    m_rootTupleSvc->addItem("tree_1", "count",  &m_count);
    m_rootTupleSvc->addItem("tree_1", "square", &m_square);
    m_rootTupleSvc->addItem("tree_1", "int",    &m_int);
    m_rootTupleSvc->addItem("tree_1", "uint",   &m_uint);
    m_rootTupleSvc->addItem("tree_1", "float",  &m_float);
    m_rootTupleSvc->addItem("tree_1", "array[2]",m_array);
    m_rootTupleSvc->addItem("tree_1", "name",    m_name);
#if 1
    // test creation of a second ROOT file
    m_rootTupleSvc->addItem("t2", "float2", &m_float2, "other.root");
#endif
    // test of a second tree in original file
    m_rootTupleSvc->addItem("tree_2","count", &m_count);
    m_rootTupleSvc->addItem("tree_2","square", &m_square);

    // test creation of memory resident tuple
    m_rootTupleSvc->addItem("memoryTree","memoryFloat",&m_memoryFloat, "", false);
    m_rootTupleSvc->addItem("memoryTree","memoryInt",&m_memoryInt,"",false);

    // check that we can find a previous item

    float* test;

    bool isFloat = m_rootTupleSvc->getItem("tree_1","float",  (void*&)test)=="Float_t";
    if( !isFloat || test!=&m_float){
        log << MSG::ERROR << "Did not retrieve a float" << endreq;
        sc = StatusCode::FAILURE;
    }else{
        log << MSG::INFO << "Found previous entry OK" << endreq;
    }

    return sc;
}


StatusCode writeJunkAlg::execute() {
     
    static int callCount = 0;

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    // note that setting these variables is all that is necessary to have it changed in the tuple itself
    ++m_count;
    m_square= m_count*m_count;
    // test int value
    m_int=(int)floor(m_count);
    m_uint=(unsigned int)floor(m_count);
    // see that array really works
    m_array[0]= m_int;
    m_array[1]= 2*m_int;

    m_float2 = m_count;

    // Set up memory resident tuple
    m_memoryFloat = m_count;
    m_memoryInt = 1;

    void *memoryTree;
    long long numMemoryEntries = m_rootTupleSvc->getOutputTreePtr(memoryTree, "memoryTree");
    log << MSG::INFO << "Memory Tree Contains " << numMemoryEntries << " entries" << endreq;

    void *testFloatPtr;
    m_rootTupleSvc->getItem("memoryTree","memoryFloat",testFloatPtr);
    float testFloat = *reinterpret_cast<float*>(testFloatPtr);

    log << MSG::INFO << "Retrieved float from memory resident tree: " << testFloat << endreq;

    // test adding names
    static const char* names[] = {"name1", "name2"};
    strncpy(m_name, (m_count<3? names[0]: names[1]) , 10);

    // check finite test routines
    m_float = m_count==5? m_count/0.0 : m_count;

    // Test the ability to turn off a row
    m_rootTupleSvc->storeRowFlag("tree_1",true) ; //callCount == 5);
    m_rootTupleSvc->storeRowFlag("memoryTree",true);
    ++callCount;


    return sc;
}


StatusCode writeJunkAlg::finalize() {
    StatusCode  sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize writeJunkAlg " << endreq;
 
    return StatusCode::SUCCESS;
}

