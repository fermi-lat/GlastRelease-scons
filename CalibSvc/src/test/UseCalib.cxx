//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibData/CalibTest1.h"
#include "CalibSvc/ICalibPathSvc.h"

/// Simple algorithm to test functioning of "the other" TDS
class UseCalib : public Algorithm {

public:
  UseCalib(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  IDataProviderSvc* m_pCalibDataSvc;
  ICalibPathSvc*    m_pCalibPathSvc;

  std::string       m_path;      // will hold path string for TDS

};

/// Instantiation of a static factory to create instances of this algorithm
//static const AlgFactory<UseCalib> Factory;
//const IAlgFactory& UseCalibFactory = Factory;
DECLARE_ALGORITHM_FACTORY(UseCalib);


UseCalib::UseCalib( const std::string&  name, 
		    ISvcLocator*        pSvcLocator )
  : Algorithm     ( name, pSvcLocator ), m_pCalibDataSvc(0)
{
  // Declare properties here.

}


StatusCode UseCalib::initialize() {
  StatusCode sc;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initialize()" << endreq;

  // So far don't have any properties, but in case we do some day..
  setProperties();


  sc = service("CalibDataSvc", m_pCalibDataSvc, true);

  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not get IDataProviderSvc interface of CalibDataSvc" 
	<< endreq;
    return sc;
  } else {
    log << MSG::DEBUG 
	<< "Retrieved IDataProviderSvc interface of CalibDataSvc" 
	<< endreq;
  }

  sc = service("CalibDataSvc", m_pCalibPathSvc, true);

  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not get ICalibPathSvc interface of CalibDataSvc" 
	<< endreq;
    return sc;
  }

  m_path = 
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_CalibTest1,
                                  std::string("vanilla") );

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;

}


StatusCode UseCalib::execute( ) {

  MsgStream log(msgSvc(), name());

  //  std::string fullPath = "/Calib/Test_1/vanilla";

  SmartDataPtr<CalibData::CalibTest1> test1Copy(m_pCalibDataSvc, m_path);

  if (!test1Copy) {
    log << MSG::ERROR << "Failed access to CalibTest1 via smart ptr" << endreq;
    return StatusCode::FAILURE;
  }

  log << MSG::INFO 
      << "Test_1 obj, serial #" <<  test1Copy->getSerNo() 
      << "  value = " << test1Copy->getValue() << endreq;
  log << MSG::INFO << "Vstart: " <<  (test1Copy->validSince()).hour(true)
      << "  Vend: " << (test1Copy->validTill()).hour(true) << endreq;

  m_pCalibDataSvc->updateObject((CalibData::CalibTest1 *)test1Copy);

  if (!test1Copy) {
    log << MSG::ERROR 
        << "Failed access to CalibTest1 via smart ptr" << endreq;
    return StatusCode::FAILURE;
  }
  log << MSG::INFO 
      << "After update Test_1 object, serial #" <<  test1Copy->getSerNo() 
      << " has value = " << test1Copy->getValue() << endreq;
  log << MSG::INFO << "Vstart: " <<  (test1Copy->validSince()).hour(true)
      << "  Vend: " << (test1Copy->validTill()).hour(true) << endreq;

  return StatusCode::SUCCESS;
}

StatusCode UseCalib::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "        UseCalib FINALIZE!! "
      << endreq;
  
  return StatusCode::SUCCESS;
}

