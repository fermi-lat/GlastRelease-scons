#define ntupleWriter_writeJunkAlg_CPP 


#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"


#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/StatusCode.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"

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
    double m_count; // special value to put into tuples
    double m_square; // another

    INTupleWriterSvc *m_rootTupleSvc;

};

static const AlgFactory<writeJunkAlg>  Factory;
const IAlgFactory& writeJunkAlgFactory = Factory;

writeJunkAlg::writeJunkAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator){

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
    m_rootTupleSvc->addItem("","count", &m_count);
    m_rootTupleSvc->addItem("","square", &m_square);

    // test of a second tree in the same file
    m_rootTupleSvc->addItem("tree_2","count", &m_count);
    m_rootTupleSvc->addItem("tree_2","square", &m_square);
    return sc;
}


StatusCode writeJunkAlg::execute() {
     
    static int callCount = 0;

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    // note that setting these variables is all that is necessary to have it changed in the tuple itself
    ++m_count;
    m_square= m_count*m_count;
    // Test the ability to turn off a row

    if (callCount != 5)  { m_rootTupleSvc->storeRowFlag(true);}
    ++callCount;


    return sc;
}


StatusCode writeJunkAlg::finalize() {
    StatusCode  sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize writeJunkAlg " << endreq;
 
    return StatusCode::SUCCESS;
}

