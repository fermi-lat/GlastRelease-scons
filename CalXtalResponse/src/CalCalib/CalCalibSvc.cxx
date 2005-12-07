// Implementation file for CalCalibSvc
//
// $Header $
//
// Author: Zachary Fewtrell

// LOCAL
#include "CalCalibSvc.h"

// GLAST

// EXTLIB
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Incident.h"

// STD

using namespace CalUtil;

static SvcFactory< CalCalibSvc > a_factory;
const ISvcFactory& CalCalibSvcFactory = a_factory; 

CalCalibSvc::CalCalibSvc(const string& name, ISvcLocator* Svc) 
  : Service(name,Svc),
    m_calibDataSvc(0),
    m_dataProviderSvc(0)
    
{
  // declare the properties
  declareProperty("CalibDataSvc",      m_calibDataSvcName = 
                  "CalibDataSvc");
  declareProperty("idealCalibXMLPath", m_idealCalibXMLPath = 
                  "$(CALXTALRESPONSEROOT)/xml/idealCalib_flight.xml");
  declareProperty("DefaultFlavor", m_defaultFlavor    
                  = "ideal");
  declareProperty("FlavorIntNonlin", m_flavorIntNonlin  = "");
  declareProperty("FlavorAsym",      m_flavorAsym       = "");
  declareProperty("FlavorPed",       m_flavorPed        = "");
  declareProperty("FlavorMeVPerDac", m_flavorMPD        = "");
  declareProperty("FlavorTholdCI",   m_flavorTholdCI    = "");
  declareProperty("FlavorTholdMuon", m_flavorTholdMuon  = "");
}

StatusCode  CalCalibSvc::queryInterface (const IID& riid, void **ppvIF) {
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
  if (!m_flavorTholdMuon.value().length())  m_flavorTholdMuon = m_defaultFlavor;

  msglog << MSG::DEBUG << "Initializing..."     << endreq;
  msglog << MSG::DEBUG << "  CalibDavaSvc   : " << m_calibDataSvcName  << endreq;
  msglog << MSG::DEBUG << "  DefaultFlavor  : " << m_defaultFlavor     << endreq;
  msglog << MSG::DEBUG << "  FlavorIntNonlin: " << m_flavorIntNonlin   << endreq;
  msglog << MSG::DEBUG << "  FlavorAsym     : " << m_flavorAsym        << endreq;
  msglog << MSG::DEBUG << "  FlavorPed      : " << m_flavorPed         << endreq;
  msglog << MSG::DEBUG << "  FlavorMeVPerDac: " << m_flavorMPD         << endreq;  
  msglog << MSG::DEBUG << "  FlavorTholdCI  : " << m_flavorTholdCI     << endreq;    
  msglog << MSG::DEBUG << "  FlavorTholdMuon: " << m_flavorTholdMuon   << endreq;   

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

  // Load ideal flavor values from xml file
  sc = loadIdealCalib();
  if (sc.isFailure()) return sc;

  // Initialize individual CalibItemMgr members.
  m_mpdMgr.initialize(m_flavorMPD,             *this);
  m_pedMgr.initialize(m_flavorPed,             *this);
  m_asymMgr.initialize(m_flavorAsym,           *this);
  m_intNonlinMgr.initialize(m_flavorIntNonlin, *this);
  m_tholdMuonMgr.initialize(m_flavorTholdMuon, *this);
  m_tholdCIMgr.initialize(m_flavorTholdCI,     *this);


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
    m_intNonlinMgr.invalidate();
    m_tholdMuonMgr.invalidate();
    m_tholdCIMgr.invalidate();
  }
  return; 
}



StatusCode CalCalibSvc::evalFaceSignal(RngIdx rngIdx, float adcPed, float &ene) {
  StatusCode sc;

  // adc -> dac
  float dac;
  sc = evalDAC(rngIdx, adcPed, dac);
  if (sc.isFailure()) return sc;

  // CalXtalID w/ no face or range info
  XtalIdx xtalIdx(rngIdx.getXtalIdx());

  // MeVPerDAC
  ValSig mpdLrg, mpdSm;
  // need to create tmp rngIdx w/ only twr/lyr/col info
  sc = getMPD(xtalIdx, mpdLrg, mpdSm);
  if (sc.isFailure()) return sc;

  RngNum rng(rngIdx.getRng());

  float mpd;
  float asymCtr;

  if (rng.getDiode() == LRG_DIODE) {
    mpd = mpdLrg.getVal();
    sc = evalAsymLrg(xtalIdx, 0, asymCtr);
    if (sc.isFailure()) return sc;
  }
  else { // diode == SM_DIODE
    mpd = mpdSm.getVal();
    sc = evalAsymSm(xtalIdx, 0, asymCtr);
    if (sc.isFailure()) return sc;
  }

  // 1st multiply each dac val by overall gain.
  ene = dac*mpd;

  // 2nd correct for overally asymmetry of diodes (use asym at center
  // of xtal)
  if (rngIdx.getFace() == POS_FACE)
    ene *= exp(-1*asymCtr/2);
  else // face == NEG_FACE
    ene *= exp(asymCtr/2);

  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::getMPD(XtalIdx xtalIdx, 
                               CalArray<DiodeNum, float> &mpd){
  ValSig mpdLrg, mpdSm;
  StatusCode sc;

  sc = getMPD(xtalIdx, mpdLrg, mpdSm);
  if (sc.isFailure()) return sc;

  mpd[LRG_DIODE] = mpdLrg.getVal();
  mpd[SM_DIODE]  = mpdSm.getVal();

  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::getPed(XtalIdx xtalIdx,
                               CalArray<XtalRng, float> &peds,
                               CalArray<XtalRng, float> &sigs) {
  float ped, sig, cos;
  StatusCode sc;

  for (XtalRng xRng; xRng.isValid(); xRng++) {
    RngIdx rngIdx(xtalIdx,
                  xRng.getFace(),
                  xRng.getRng());

    sc = getPed(rngIdx, ped, sig, cos);
    if (sc.isFailure()) return sc;

    peds[xRng] = ped;
    sigs[xRng] = sig;
  }

  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::getTholdCI(XtalIdx xtalIdx,
                                   CalArray<XtalDiode, float> &trigThresh,
                                   CalArray<FaceNum, float> &lacThresh){
  StatusCode sc;
  ValSig fle, fhe, lac;

  for (FaceNum face; face.isValid(); face++) {
    FaceIdx faceIdx(xtalIdx,face);
    
    sc = getTholdCI(faceIdx, fle, fhe, lac);
    if (sc.isFailure()) return sc;

    trigThresh[XtalDiode(face, LRG_DIODE)] = fle.getVal();
    trigThresh[XtalDiode(face, SM_DIODE)]  = fhe.getVal();
    lacThresh[face] = lac.getVal();
  }
  
  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::getULDCI(XtalIdx xtalIdx,
                                 CalArray<XtalRng, float> &uldThold){
  StatusCode sc;
  ValSig uld;

  for (XtalRng xRng; xRng.isValid(); xRng++) {
    RngIdx rngIdx(xtalIdx, xRng);

    sc = getULDCI(rngIdx, uld);
    if (sc.isFailure()) return sc;

    uldThold[xRng] = uld.getVal();
  }

  return StatusCode::SUCCESS;
}
