//$Header$
#include <stdio.h>

// #include "CalibData/CalibTime.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
// #include "GaudiKernel/IDetDataSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibData/CalibTest1.h"
#include "CalibData/CalibModel.h"


/// Simple algorithm to test functioning of "the other" TDS
class UseCalib : public Algorithm {

public:
  UseCalib(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  IDataProviderSvc* m_pCalibDataSvc;
  // Maybe something to say which kind of data to look up?

};

/// Instantiation of a static factory to create instances of this algorithm
static const AlgFactory<UseCalib> Factory;
const IAlgFactory& UseCalibFactory = Factory;


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


  //  IDataProviderSvc* m_pCalibDataSvc;
  sc = service("CalibDataSvc", m_pCalibDataSvc, true);

  // Query the IDataProviderSvc interface of the calib data service
  /*  sc = calibSvc->queryInterface(IID_IDataProviderSvc, 
                                (void**) &m_pCalibDataSvc);
  */
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

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  /*
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }
  log << MSG::DEBUG << "Properties were read from jobOptions" << endreq;
  log << MSG::INFO << "Time of first event: (seconds since 1970) "
      << m_startTime 
      << std::dec << ")" 
      << endreq; 
  log << MSG::INFO << "Time between events (seconds): "
      << m_delayTime 
      << endreq;

  */
  return StatusCode::SUCCESS;

}


StatusCode UseCalib::execute( ) {

  MsgStream log(msgSvc(), name());

  //  SmartDataPtr<CalibData::CalibTest1> test1(m_pCalibDataSvc,
  //                                        CalibData::Test_Gen);
  //  CalibData::CalibTest1* test1 = 
  //    SmartDataPtr<CalibData::CalibTest1>(m_pCalibDataSvc, CalibData::Test_Gen);
  
  
  std::string fullPath = CalibData::Test_1 + "/vanilla";
  DataObject *pObject;

  m_pCalibDataSvc->retrieveObject(fullPath, pObject);
  //  m_pCalibDataSvc->retrieveObject(CalibData::Test_1, test1Copy);

  CalibData::CalibTest1* test1Copy = 0;
  try {
    test1Copy = dynamic_cast<CalibData::CalibTest1 *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR << "Dynamic cast to CalibTest1 failed" << endreq;
    return StatusCode::FAILURE;
  }
  log << MSG::INFO 
      << "Test_1 object, serial #" <<  test1Copy->getSerNo() 
      << " has value name of "  << test1Copy->getValueName() 
      << " and value = " << test1Copy->getValue() << endreq;

  return StatusCode::SUCCESS;
}



StatusCode UseCalib::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "------------- FINALIZE!! -------------------------------------------"
      << endreq;
  
  return StatusCode::SUCCESS;
}

