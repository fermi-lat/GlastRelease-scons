// $Header$

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "CalUtil/CalFailureModeSvc.h"
#include "idents/CalXtalId.h"

// Define the class here instead of in a header file: not needed anywhere but here!
//------------------------------------------------------------------------------
/** 
A simple algorithm.

  
*/
class test_CalUtil : public Algorithm {
public:
    test_CalUtil(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 

  //! number of times called
  int m_count; 

  /// pointer to failure mode service
  ICalFailureModeSvc* m_FailSvc;

};
//------------------------------------------------------------------------

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( test_CalUtil );

static const AlgFactory<test_CalUtil>  Factory;
const IAlgFactory& test_CalUtilFactory = Factory;

//------------------------------------------------------------------------
//! ctor
test_CalUtil::test_CalUtil(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0)
{
}

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode test_CalUtil::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    sc = service("CalFailureModeSvc", m_FailSvc);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to find CalFailureMode service" << endreq;
        return sc;
    }

    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode test_CalUtil::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "executing " << ++m_count << " time" << endreq;

    idents::CalXtalId id1(10,2,2);  // tower 10, layer 3(y)
    idents::CalXtalId id2(11,1,2);  // tower 11, layer 1(y)
    idents::CalXtalId id3(3,5,3);   // tower 3, layer 5(y)
    idents::CalXtalId id4(4,6,3);   // tower 4, layer 6(x)
    idents::CalXtalId id5(4,1,3);   // tower 4, layer 1(y)

    if (m_FailSvc == 0) return StatusCode::FAILURE;
    if (m_FailSvc->matchChannel(id1,idents::CalXtalId::NEG)) {
      log << MSG::INFO << "removed channel (10,2,2) NEG" << endreq;
    } else {
      log << MSG::ERROR << "failed to remove channel (10,2,2) NEG" << endreq;
      return StatusCode::FAILURE;
    }
    if (m_FailSvc->matchChannel(id2,idents::CalXtalId::NEG)) {
      log << MSG::INFO << "removed channel (11,1,2) NEG" << endreq;
    } else {
      log << MSG::ERROR << "failed to remove channel (11,1,2) NEG" << endreq;
      return StatusCode::FAILURE;
    }
    if (!m_FailSvc->matchChannel(id3,idents::CalXtalId::NEG)) {
      log << MSG::INFO << "left channel (3,5,3) NEG" << endreq;
    } else {
      log << MSG::ERROR << "erroneously removed channel (3,5,3) NEG" << endreq;
      return StatusCode::FAILURE;
    }
    if (m_FailSvc->matchChannel(id4,idents::CalXtalId::NEG)) {
      log << MSG::INFO << "removed channel (4,6,3) NEG" << endreq;
    } else {
      log << MSG::ERROR << "erroneously left channel (4,6,3) NEG" << endreq;
      return StatusCode::FAILURE;
    }
    if (!m_FailSvc->matchChannel(id5,idents::CalXtalId::NEG)) {
      log << MSG::INFO << "left channel (4,1,3) NEG" << endreq;
    } else {
      log << MSG::ERROR << "erroneously removed channel (4,1,3) NEG" << endreq;
      return StatusCode::FAILURE;
    }


    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode test_CalUtil::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize after " << m_count << " calls." << endreq;
    
    return sc;
}