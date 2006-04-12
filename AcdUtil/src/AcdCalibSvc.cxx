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

// EXTLIB
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Incident.h"

// STD

static SvcFactory< AcdCalibSvc > a_factory;
const ISvcFactory& AcdCalibSvcFactory = a_factory; 

AcdCalibSvc::AcdCalibSvc(const std::string& name, ISvcLocator* Svc) 
  : Service(name,Svc),
    m_calibDataSvc(0),
    m_dataProviderSvc(0)    
{
  // declare the properties
  declareProperty("CalibDataSvc",    m_calibDataSvcName = "CalibDataSvc");
  declareProperty("DefaultFlavor",   m_defaultFlavor    = "ideal");
  declareProperty("FlavorPed",       m_flavorPed        = "");
  declareProperty("FlavorMipPeak",   m_flavorMipPeak    = "");
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

  // Post-init processing of jobOptions.
  // set default flavors unless otherwise specified
  if (!m_flavorPed.value().length())        m_flavorPed       = m_defaultFlavor;
  if (!m_flavorMipPeak.value().length())    m_flavorMipPeak   = m_defaultFlavor;

  msglog << MSG::DEBUG << "Initializing..."     << endreq;
  msglog << MSG::DEBUG << "  CalibDavaSvc   : " << m_calibDataSvcName  << endreq;
  msglog << MSG::DEBUG << "  DefaultFlavor  : " << m_defaultFlavor     << endreq;
  msglog << MSG::DEBUG << "  FlavorPed      : " << m_flavorPed         << endreq;
  msglog << MSG::DEBUG << "  FlavorMipPeak  : " << m_flavorMipPeak     << endreq;  

  // Grab pointer to CalibDataSvc
  sc = service(m_calibDataSvcName, m_calibDataSvc, true);
  if ( !sc.isSuccess() ) {
    msglog << MSG::ERROR << "Could not get CalibDataSvc" << endreq;
    return sc;
  }

  // Query the IDataProvider interface of the CalibDataService
  sc = m_calibDataSvc->queryInterface(IID_IDataProviderSvc, 
                                      (void**) &m_dataProviderSvc);
  if ( !sc.isSuccess() ) {
    msglog << MSG::ERROR 
           << "Could not query IDataProviderSvc interface of CalibDataSvc" 
           << endreq;
    return sc;
  }

  // Initialize individual CalibItemMgr members.
  m_pedMgr.initialize(m_flavorPed,             *this);
  m_gainMgr.initialize(m_flavorMipPeak,        *this);

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

/// Inform that a new incident has occured
void AcdCalibSvc::handle ( const Incident& inc ) { 
  if ((inc.type() == "BeginEvent")) {
    m_pedMgr.invalidate();
    m_gainMgr.invalidate();
  }
  return; 
}

