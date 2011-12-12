// $Header$
/** @file
    @author Z.Fewtrell
*/

// LOCAL INCLUDES
#include "CalCalibShared.h"

// GLAST INCLUDES

// EXTLIB INCLUDES
#include "GaudiKernel/MsgStream.h"

// STD INCLUDES

/// intialize / retrieve all Gaudi based member values
StatusCode CalCalibShared::initialize(Service &service) {
  StatusCode sc;

  m_service = &service;

  sc = m_idealCalib.readCfgFile(m_idealCalibXMLPath);
  if (sc.isFailure()) return sc;

  IService *calibDataSvc;
  // Grab pointer to CalibDataSvc
  sc = m_service->service(m_calibDataSvcName, calibDataSvc, true);
  if ( !sc.isSuccess() ) {
    MsgStream msglog(m_service->msgSvc(), m_service->name()); 
    msglog << MSG::ERROR << "Could not get CalibDataSvc" << endreq;
    return sc;
  }

  // Query the IDataProvider interface of the CalibDataService
  sc = calibDataSvc->queryInterface(IDataProviderSvc::interfaceID(), 
                                    (void**) &m_dataProviderSvc);
  if ( !sc.isSuccess() ) {
    MsgStream msglog(m_service->msgSvc(), m_service->name()); 
    msglog << MSG::ERROR 
           << "Could not query IDataProviderSvc interface of CalibDataSvc" 
           << endreq;
    return sc;
  }

  // grab pointer to DataPathSvc
  sc = m_service->service(m_calibDataSvcName, m_calibPathSvc, true);
  if ( !sc.isSuccess() ) {
    MsgStream msglog(m_service->msgSvc(), m_service->name());
    msglog << MSG::ERROR << "Could not get CalibDataSvc" << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
  
}
