#define ntupleWriter_writeJunkAlg_CPP 

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"


#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/INTuple.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/StatusCode.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"

//------------------------------------------------------------------------------
/*! \class writeJunkAlg
\brief test algorithm for the ntupleWriterSvc

This algorithm tests the creation and writing of ntuples via
the ntupleWriterSvc.  The output from this routine is a
ROOT ntuple, containing 19 entries.
*/

class writeJunkAlg : public Algorithm {
    
public:
    //! Constructor of this form must be provided
    writeJunkAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    std::string m_tupleName;
    float m_float;

    INTupleWriterSvc *m_ntupleWriteSvc;

};

//------------------------------------------------------------------------------
static const AlgFactory<writeJunkAlg>  Factory;
const IAlgFactory& writeJunkAlgFactory = Factory;
//------------------------------------------------------------------------------
/// 
writeJunkAlg::writeJunkAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator){

    declareProperty("tupleName",  m_tupleName="");

}

//------------------------------------------------------------------------------
/*! 
*/
StatusCode writeJunkAlg::initialize() {
    
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    // get a pointer to our ntupleWriterSvc
    sc = service("ntupleWriterSvc", m_ntupleWriteSvc);

    if( sc.isFailure() ) {
        log << MSG::ERROR << "writeJunkAlg failed to get the ntupleWriterSvc" << endreq;
        return sc;
    }

    return sc;
}



//------------------------------------------------------------------------------
StatusCode writeJunkAlg::execute() {
     
    static int callCount = 0;

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "CallCount", callCount);

    m_float = 700.05f;
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "MyFirstItem", m_float);

    float zero = 0.0;
    // Let's try to add an undefined float
    float bad = 100.0/zero;
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "BadValue", bad);

    // Test the ability to turn off a row
    if (callCount == 5) m_ntupleWriteSvc->storeRowFlag(false);
    // test call to store ntuples to disk during execution
    if (callCount == 6) m_ntupleWriteSvc->saveNTuples();

    ++callCount;

    return sc;
}


//------------------------------------------------------------------------------
StatusCode writeJunkAlg::finalize() {
    StatusCode  sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize writeJunkAlg " << endreq;
 
    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------



