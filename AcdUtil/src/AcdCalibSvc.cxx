// $Header $
/** @file
    @author Zach Fewtrell
 */
// @file
//
//
// Author: Zachary Fewtrell

// LOCAL
#include "AcdCalibSvc.h"

// GLAST
#include "CalibData/Acd/AcdCalib.h"

// EXTLIB
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Incident.h"

// STD

//static SvcFactory< AcdCalibSvc > a_factory;
//const ISvcFactory& AcdCalibSvcFactory = a_factory; 
DECLARE_SERVICE_FACTORY(AcdCalibSvc);

AcdCalibSvc::AcdCalibSvc(const std::string& name, ISvcLocator* Svc) 
  : Service(name,Svc),
    m_calibDataSvc(0),
    m_dataProviderSvc(0)    
{
  // declare the properties
  declareProperty("CalibDataSvc",    m_calibDataSvcName = "CalibDataSvc");
  declareProperty("DefaultFlavor",   m_defaultFlavor    = "ideal");
  addCalibration( new AcdCalibMgrTmpl< CalibData::AcdPedCalib >,"FlavorPed");
  addCalibration( new AcdCalibMgrTmpl< CalibData::AcdGainCalib >,"FlavorGain");
  addCalibration( new AcdCalibMgrTmpl< CalibData::AcdHighRangeCalib >,"FlavorHighRange");
  addCalibration( new AcdCalibMgrTmpl< CalibData::AcdCoherentNoiseCalib >,"FlavorCoherentNoise");
}

StatusCode  AcdCalibSvc::queryInterface (const InterfaceID& riid, void **ppvIF) {
  if (IID_IAcdCalibSvc == riid) {
    *ppvIF = (IAcdCalibSvc*)(this);
    return StatusCode::SUCCESS;
  } else return Service::queryInterface (riid, ppvIF);
}

StatusCode AcdCalibSvc::initialize () 
{
  // Call super-class
  Service::initialize ();

  MsgStream msglog(msgSvc(), name()); 
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  //-- jobOptions --//
  if ((sc = setProperties()).isFailure()) {
    msglog << MSG::ERROR << "Failed to set properties" << endreq;
    return sc;
  }

  msglog << MSG::DEBUG << "Initializing..."     << endreq;
  msglog << MSG::DEBUG << "  CalibDavaSvc         : " << m_calibDataSvcName  << endreq;
  msglog << MSG::DEBUG << "  DefaultFlavor        : " << m_defaultFlavor     << endreq;

  // Grab pointer to CalibDataSvc
  sc = service(m_calibDataSvcName, m_calibDataSvc, true);
  if ( !sc.isSuccess() ) {
    msglog << MSG::ERROR << "Could not get CalibDataSvc" << endreq;
    return sc;
  }
  
  sc = service("CalibDataSvc", m_calibPathSvc, true);
  if ( !sc.isSuccess() ) {
    msglog << MSG::ERROR << "Could not get CalibPathSvc" << endreq;
    return sc;
  }

  // Query the IDataProvider interface of the CalibDataService
  sc = m_calibDataSvc->queryInterface(IDataProviderSvc::interfaceID(), 
                                      (void**) &m_dataProviderSvc);
  if ( !sc.isSuccess() ) {
    msglog << MSG::ERROR 
           << "Could not query IDataProviderSvc interface of CalibDataSvc" 
           << endreq;
    return sc;
  }

  sc = prepapreManagers(msglog,m_defaultFlavor);
  if ( !sc.isSuccess() ) {
    msglog << MSG::ERROR 
           << "Could not prepare calibration managers"
           << endreq;
    return sc;
  }

  // Get ready to listen for BeginEvent
  IIncidentSvc* incSvc;
  sc = service("IncidentSvc", incSvc, true);
  if (sc.isSuccess() ) {
    int priority = 50;  // this should be lower priority (higher #?) than CalibDataSvc
    incSvc->addListener(this, "BeginEvent", priority);
  } else {
    msglog << MSG::ERROR << "can't find IncidentSvc" << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}

const std::string AcdCalibSvc::getCalibPath(const ICalibPathSvc::CalibItem item, const std::string& flavor) const {
  return m_calibPathSvc->getCalibPath(item, flavor);
}

void AcdCalibSvc::addCalibration(AcdCalibMgr* calibMgr,const std::string& flavorName) {
  AcdCalibData::CALTYPE calibType = calibMgr->calibType();
  StringProperty* property = new StringProperty(flavorName,"");
  declareProperty(flavorName,*property);
  addMgr(calibType, calibMgr, property);
}
