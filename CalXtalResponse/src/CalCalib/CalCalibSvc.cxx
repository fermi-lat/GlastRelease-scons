// $Header $
/** @file
    @author Zach Fewtrell
 */
// @file
//
//
// Author: Zachary Fewtrell

// LOCAL
#include "CalCalibSvc.h"

// GLAST

// EXTLIB
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"

// STD

using namespace CalUtil;

static SvcFactory< CalCalibSvc > a_factory;
const ISvcFactory& CalCalibSvcFactory = a_factory; 

CalCalibSvc::CalCalibSvc(const string& name, ISvcLocator* Svc) 
  : Service(name,Svc),
    m_pedMgr(m_ccsShared),
    m_inlMgr(m_ccsShared),
    m_asymMgr(m_ccsShared),
    m_mpdMgr(m_ccsShared),
    m_tholdCIMgr(m_ccsShared)
{
  // declare the properties
  declareProperty("CalibDataSvc",      m_ccsShared.m_calibDataSvcName = 
                  "CalibDataSvc");
  declareProperty("idealCalibXMLPath", m_ccsShared.m_idealCalibXMLPath = 
                  "$(CALXTALRESPONSEROOT)/xml/idealCalib_flight.xml");
  declareProperty("DefaultFlavor", m_defaultFlavor    
                  = "ideal");
  declareProperty("FlavorIntNonlin", m_flavorIntNonlin  = "");
  declareProperty("FlavorAsym",      m_flavorAsym       = "");
  declareProperty("FlavorPed",       m_flavorPed        = "");
  declareProperty("FlavorMeVPerDac", m_flavorMPD        = "");
  declareProperty("FlavorTholdCI",   m_flavorTholdCI    = "");
}

StatusCode  CalCalibSvc::queryInterface (const InterfaceID& riid, void **ppvIF) {
  if (IID_ICalCalibSvc == riid) {
    *ppvIF = (ICalCalibSvc*)(this);
    return StatusCode::SUCCESS;
  } else return Service::queryInterface (riid, ppvIF);
}

StatusCode CalCalibSvc::initialize () 
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
  if (!m_flavorIntNonlin.value().length())  m_flavorIntNonlin = m_defaultFlavor;
  if (!m_flavorAsym.value().length())       m_flavorAsym      = m_defaultFlavor;
  if (!m_flavorPed.value().length())        m_flavorPed       = m_defaultFlavor;
  if (!m_flavorMPD.value().length())        m_flavorMPD       = m_defaultFlavor;
  if (!m_flavorTholdCI.value().length())    m_flavorTholdCI   = m_defaultFlavor;

  msglog << MSG::DEBUG << "Initializing..."     << endreq;
  msglog << MSG::DEBUG << "  CalibDavaSvc   : " << m_ccsShared.m_calibDataSvcName  << endreq;
  msglog << MSG::DEBUG << "  DefaultFlavor  : " << m_defaultFlavor     << endreq;
  msglog << MSG::DEBUG << "  FlavorIntNonlin: " << m_flavorIntNonlin   << endreq;
  msglog << MSG::DEBUG << "  FlavorAsym     : " << m_flavorAsym        << endreq;
  msglog << MSG::DEBUG << "  FlavorPed      : " << m_flavorPed         << endreq;
  msglog << MSG::DEBUG << "  FlavorMeVPerDac: " << m_flavorMPD         << endreq;  
  msglog << MSG::DEBUG << "  FlavorTholdCI  : " << m_flavorTholdCI     << endreq;    


  /// init all the CalCalibShared stuff
  sc = m_ccsShared.initialize(*this);
  if (sc.isFailure()) return sc;

  // Initialize individual CalibItemMgr members.
  sc = m_mpdMgr.initialize(m_flavorMPD);
  if (sc.isFailure()) return sc;
  sc = m_pedMgr.initialize(m_flavorPed);
  if (sc.isFailure()) return sc;
  sc = m_asymMgr.initialize(m_flavorAsym);
  if (sc.isFailure()) return sc;
  sc = m_inlMgr.initialize(m_flavorIntNonlin);
  if (sc.isFailure()) return sc;
  sc = m_tholdCIMgr.initialize(m_flavorTholdCI);
  if (sc.isFailure()) return sc;

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
void CalCalibSvc::handle ( const Incident& inc ) { 
  if ((inc.type() == "BeginEvent")) {
    m_mpdMgr.invalidate();
    m_pedMgr.invalidate();
    m_asymMgr.invalidate();
    m_inlMgr.invalidate();
    m_tholdCIMgr.invalidate();
  }
  return; 
}

StatusCode CalCalibSvc::evalFaceSignal(RngIdx rngIdx, float adcPed, float &ene) {
  StatusCode sc;

  // adc -> cidac
  float cidac;
  sc = evalCIDAC(rngIdx, adcPed, cidac);
  if (sc.isFailure()) return sc;

  // CalXtalID w/ no face or range info
  XtalIdx xtalIdx(rngIdx.getXtalIdx());

  // MeVPerDAC
  const CalMevPerDac *calMPD;
  // need to create tmp rngIdx w/ only twr/lyr/col info
  calMPD = getMPD(xtalIdx);
  if (!calMPD) return StatusCode::FAILURE;

  RngNum rng(rngIdx.getRng());

  float mpd;
  float asymCtr;

  if (rng.getDiode() == LRG_DIODE) {
    mpd = calMPD->getBig()->getVal();
    sc = getAsymCtr(xtalIdx, ASYM_LL, asymCtr);
    if (sc.isFailure()) return sc;
  }
  else { // diode == SM_DIODE
    mpd = calMPD->getSmall()->getVal();
    sc = getAsymCtr(xtalIdx, ASYM_SS, asymCtr);
    if (sc.isFailure()) return sc;
  }

  // 1st multiply each cidac val by overall gain.
  ene = cidac*mpd;

  // 2nd correct for overally asymmetry of diodes (use asym at center
  // of xtal)
  if (rngIdx.getFace() == POS_FACE)
    ene *= exp(-1*asymCtr/2);
  else // face == NEG_FACE
    ene *= exp(asymCtr/2);

  return StatusCode::SUCCESS;
}



