//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibData/Nas/CalibLATAlignment.h"
#include "CalibSvc/ICalibPathSvc.h"

/// Simple algorithm to test functioning of "the other" TDS
class UseLATAlignment : public Algorithm {

public:
  UseLATAlignment(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  IDataProviderSvc* m_pCalibDataSvc;
  ICalibPathSvc*    m_pCalibPathSvc;
  std::string       m_path;
  // Maybe something to say which kind of data to look up?

};

/// Instantiation of a static factory to create instances of this algorithm
static const AlgFactory<UseLATAlignment> Factory;
const IAlgFactory& UseLATAlignmentFactory = Factory;


UseLATAlignment::UseLATAlignment( const std::string&  name, 
		    ISvcLocator*        pSvcLocator )
  : Algorithm     ( name, pSvcLocator ), m_pCalibDataSvc(0)
{
  // Declare properties here.

}


StatusCode UseLATAlignment::initialize() {
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
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_NAS_LATAlignment,
                                  std::string("vanilla") );

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;

}


StatusCode UseLATAlignment::execute( ) {

  MsgStream log(msgSvc(), name());

  SmartDataPtr<CalibData::CalibLATAlignment> test1Copy(m_pCalibDataSvc, m_path);

  if (!test1Copy) {
    log << MSG::ERROR << "Failed access to CalibLATAlignment via smart ptr" << endreq;
    return StatusCode::FAILURE;
  }

  double roll,pitch,yaw;  
  
  if(!(roll=test1Copy->getRoll())){
    log << MSG::ERROR 
        << "Failed to read phi" << endreq;
    return StatusCode::FAILURE;
  }  
  if(!(pitch=test1Copy->getPitch())){
    log << MSG::ERROR 
        << "Failed to read theta" << endreq;
    return StatusCode::FAILURE;
  }
  if(!(yaw=test1Copy->getYaw())){
    log << MSG::ERROR 
        << "Failed to read psi" << endreq;
    return StatusCode::FAILURE;
  }
  
  
  log << MSG::INFO 
      << "SAA boundary obj, serial #" <<  test1Copy->getSerNo() << endreq;

  log << MSG::INFO << "Vstart: " <<  (test1Copy->validSince()).hours()
      << "  Vend: " << (test1Copy->validTill()).hours() << endreq;

  log << MSG::INFO << "Roll: "   << roll  << endreq;
  log << MSG::INFO << "Pitch: " << pitch   << endreq;
  log << MSG::INFO << "Yaw: "   << yaw   << endreq;



  return StatusCode::SUCCESS;
}

StatusCode UseLATAlignment::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "        UseLATAlignment FINALIZE!! "
      << endreq;
  
  return StatusCode::SUCCESS;
}

