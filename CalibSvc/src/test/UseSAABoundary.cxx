//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibData/Nas/CalibSAABoundary.h"
#include "CalibSvc/ICalibPathSvc.h"

/// Simple algorithm to test functioning of "the other" TDS
class UseSAABoundary : public Algorithm {

public:
  UseSAABoundary(const std::string& name, ISvcLocator* pSvcLocator); 

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
//static const AlgFactory<UseSAABoundary> Factory;
//const IAlgFactory& UseSAABoundaryFactory = Factory;
DECLARE_ALGORITHM_FACTORY(UseSAABoundary);


UseSAABoundary::UseSAABoundary( const std::string&  name, 
		    ISvcLocator*        pSvcLocator )
  : Algorithm     ( name, pSvcLocator ), m_pCalibDataSvc(0)
{
  // Declare properties here.

}


StatusCode UseSAABoundary::initialize() {
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
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_NAS_SAABoundary,
                                  std::string("vanilla") );

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;

}


StatusCode UseSAABoundary::execute( ) {

  MsgStream log(msgSvc(), name());


  SmartDataPtr<CalibData::CalibSAABoundary> test1Copy(m_pCalibDataSvc, m_path);

  if (!test1Copy) {
    log << MSG::ERROR << "Failed access to CalibSAABoundary via smart ptr" << endreq;
    return StatusCode::FAILURE;
  }

  std::pair<double,double>vertex;
  
  if(!test1Copy->getFirstVertex(vertex)){
    log << MSG::ERROR 
        << "Failed to read first vertex" << endreq;
    return StatusCode::FAILURE;
  }
  
  log << MSG::INFO 
      << "SAA boundary obj, serial #" <<  test1Copy->getSerNo() << endreq;

  log << MSG::INFO << "Vstart: " <<  (test1Copy->validSince()).hour(true)
      << "  Vend: " << (test1Copy->validTill()).hour(true) << endreq;

  log << MSG::INFO << "Vertex: " <<  vertex.first
      << " , " << vertex.second << endreq;

  while(test1Copy->getNextVertex(vertex)){
      log << MSG::INFO << "Vertex: " <<  vertex.first
	  << " , " << vertex.second << endreq;        
  };
  
/*  
  m_pCalibDataSvc->updateObject((CalibData::CalibSAABoundary *)test1Copy);

  if (!test1Copy) {
    log << MSG::ERROR 
        << "Failed access to CalibTest1 via smart ptr" << endreq;
    return StatusCode::FAILURE;
  }
  log << MSG::INFO 
      << "After update Test_1 object, serial #" <<  test1Copy->getSerNo() 
      << " has value = " << test1Copy->getValue() << endreq;
  log << MSG::INFO << "Vstart: " <<  (test1Copy->validSince()).hours()
      << "  Vend: " << (test1Copy->validTill()).hours() << endreq;
*/
  return StatusCode::SUCCESS;
}

StatusCode UseSAABoundary::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "        UseSAABoundary FINALIZE!! "
      << endreq;
  
  return StatusCode::SUCCESS;
}

