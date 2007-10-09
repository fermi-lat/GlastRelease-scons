// $Header $
/** @file
    @author Zach Fewtrell
*/
// @file
//
//
// Author: Zachary Fewtrell

// LOCAL
#include "CalSignalTool.h"
#include "IXtalSignalTool.h"
#include "../Calib/IPrecalcCalibTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"

// GLAST
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"
#include "CalUtil/CalGeom.h"

// EXTLIB
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"

#include "CLHEP/Random/RandGauss.h"


// STD
#include <algorithm>
#include <utility>
#include <map>
#include <string>

static ToolFactory< CalSignalTool > a_factory;
const IToolFactory& CalSignalToolFactory = a_factory; 

CalSignalTool::CalSignalTool(const std::string& type,
                             const std::string& name,
                             const IInterface* parent) 
  : AlgTool(type,name,parent),
    m_calSignalMap(CalUtil::DiodeIdx::N_VALS),
    m_isValid(false),
    m_evtSvc(0),
    m_eTowerCAL(0),
    m_eLATTowers(0),
    m_eMeasureX(0),
    m_eXtal(0),
    m_xtalSignalTool(0),
    m_precalcCalib(0),
    m_calCalibSvc(0),
    m_detSvc(0)
{
  declareInterface<ICalSignalTool>(this);

  declareProperty("CalCalibSvc",         m_calCalibSvcName    = "CalCalibSvc");
  declareProperty("XtalSignalToolName",  m_xtalSignalToolName = "XtalSignalTool");
  declareProperty("enableNoise",         m_enableNoise = true);
  declareProperty("PrecalcCalibTool",    m_precalcCalibName   = "PrecalcCalibTool");

}

StatusCode CalSignalTool::initialize() {
  MsgStream msglog(msgSvc(), name());

  msglog << MSG::INFO << "CalSignalTool is initializing" << endreq;

  //-- jobOptions --//
  StatusCode sc;
  if ((sc = setProperties()).isFailure()) {
    msglog << MSG::ERROR << "Failed to set properties" << endreq;
    return sc;
  }

  //-- Retreive EventDataSvc
  sc = serviceLocator()->service( "EventDataSvc", m_evtSvc, true );
  if(sc.isFailure()){
    msglog << MSG::ERROR << "Could not find EventDataSvc" << endreq;
    return sc;
  }

  // now try to find the GlastDetSvc service
  sc = service("GlastDetSvc", m_detSvc);
  if (sc.isFailure() ) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::ERROR << "  Unable to get GlastDetSvc " << endreq;
    return sc;
  }

  //-- Retrieve constants from GlastDetSVc
  sc = retrieveConstants();
  if (sc.isFailure())
    return sc;

  //-- find out which tems are installed.
  m_twrList = CalUtil::findActiveTowers(*m_detSvc);

  //-- Retrieve xtalSignalTool
  sc = toolSvc()->retrieveTool("XtalSignalTool", 
                               m_xtalSignalToolName, 
                               m_xtalSignalTool, 
                               this ); 
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to create " << m_xtalSignalToolName << endreq;
    return sc;
  }

  // Get ready to listen for BeginEvent
  IIncidentSvc* incSvc;
  sc = service("IncidentSvc", incSvc, true);
  if (sc.isSuccess() ) {
    incSvc->addListener(this, "BeginEvent");
  } else {
    msglog << MSG::ERROR << "can't find IncidentSvc" << endreq;
    return sc;
  }

  // this tool may also be shared by CalTrigTool, global ownership
  sc = toolSvc()->retrieveTool("PrecalcCalibTool", 
                               m_precalcCalibName, 
                               m_precalcCalib,
                               0); // shared by other code
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to create " << m_precalcCalibName << endreq;
    return sc;
  }

  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get " << m_calCalibSvcName  << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}


StatusCode CalSignalTool::retrieveConstants() {
  double value;
  typedef map<int*,string> PARAMAP;


  PARAMAP param;
  param[&m_eTowerCAL]  =    string("eTowerCAL");
  param[&m_eLATTowers] =    string("eLATTowers");
  param[&m_eXtal]      =    string("eXtal");
  param[&m_eMeasureX]  =    string("eMeasureX");
  param[m_ePerMeV+1]     = string("cal.ePerMeVSmall");
  param[m_ePerMeV]       = string("cal.ePerMevLarge");



  for(PARAMAP::iterator it=param.begin(); it!=param.end();it++){
    if(!m_detSvc->getNumericConstByName((*it).second, &value)) {
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::ERROR << " constant " <<(*it).second 
             <<" not defined" << endreq;
      return StatusCode::FAILURE;
    } else *((*it).first)=(int)value;
  }

  return StatusCode::SUCCESS;
}


void CalSignalTool::newEvent() {

  fill(m_calSignalMap.begin(),
       m_calSignalMap.end(),
       0);

  m_calRelMap.clear();

  m_isValid = false;
}


const ICalSignalTool::CalSignalMap *CalSignalTool::getCalSignalMap() {
  /// check that all internal data is up2date
  StatusCode sc = syncData();
  if (sc.isFailure())
    return 0;

  return &m_calSignalMap;

}

const ICalSignalTool::CalRelationMap *CalSignalTool::getCalRelationMap() {
  /// check that all internal data is up2date
  StatusCode sc = syncData();
  if (sc.isFailure())
    return 0;

  return &m_calRelMap;
}

/// apply McIntegratingHits to diode signal levels and simulate noise
StatusCode CalSignalTool::syncData() {
  /// return early if we've already calcuated everything for this event
  if (m_isValid)
    return StatusCode::SUCCESS;

  /// sum individual mc hits
  StatusCode sc = loadSignalMaps();
  if (sc.isFailure())
    return sc;
  
  /// calc electronic noise (independent of hits)
  if (m_enableNoise) {
    sc = calcNoise();
    if (sc.isFailure())
      return sc;
  }

  /// indicate that event data is up-to-date for now
  m_isValid = true;
  
  return StatusCode::SUCCESS;
}

/// apply All Cal McIntegratingHits to crystal diode signal levels.
StatusCode CalSignalTool::loadSignalMaps() {
  // get McIntegratingHit collection. Abort if empty.
  SmartDataPtr<Event::McIntegratingHitVector> 
    McCalHits(m_evtSvc, EventModel::MC::McIntegratingHitCol );
  
  if (McCalHits == 0) {
    // create msglog only when needed for speed.
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::DEBUG; 
    if (msglog.isActive()){ 
      msglog.stream() << "no cal hits found" ;} 
    msglog << endreq;
    return StatusCode::SUCCESS;
  }

  // loop over hits - pick out CAL hits
  for (Event::McIntegratingHitVector::const_iterator it = McCalHits->begin(); 
       it != McCalHits->end(); 
       it++) {

    // get volumeId from hit
    const idents::VolumeIdentifier volId = 
      ((idents::VolumeIdentifier)(*it)->volumeID());

    //   extracting parameters from volume Id identifying as in CAL
    if ((int)volId[CalUtil::fLATObjects]   == m_eLATTowers &&
        (int)volId[CalUtil::fTowerObjects] == m_eTowerCAL) { 

      // apply hit signal to relavant diodes
      StatusCode sc = sumHit(**it);
      if (sc.isFailure())
        return sc;

      // store MCHit <> xtal relation
      sc = registerHitRel(**it);
      if (sc.isFailure())
        return sc;

    }
  }


  return StatusCode::SUCCESS;
}

/// apply single crystal hit to signal level
StatusCode CalSignalTool::sumHit(const Event::McIntegratingHit &hit) {
  // retrieve crystal id from hit
  using namespace idents;
  using namespace CalUtil;

  const XtalIdx xtalIdx(CalXtalId(hit.volumeID()));

  // destination for signal data from this hit
  CalUtil::CalArray<CalUtil::XtalDiode, float> hitSignal;
  fill(hitSignal.begin(), hitSignal.end(), 0);

  // convert the hit into signal level
  StatusCode sc = m_xtalSignalTool->calculate(hit, 
                                              hitSignal);
  if (sc.isFailure())
    return sc;

  // sum xtal signals into full arry
  for (XtalDiode xtalDiode;
       xtalDiode.isValid();
       xtalDiode++) {
    const DiodeIdx diodeIdx(xtalIdx, xtalDiode);
    
    m_calSignalMap[diodeIdx] += hitSignal[xtalDiode];
  }
  
  return StatusCode::SUCCESS;
}

/// store crystal <> McHit relation
StatusCode CalSignalTool::registerHitRel(Event::McIntegratingHit &hit) {
  /// retrieve crystal id from hit
  const CalUtil::XtalIdx xtalIdx(idents::CalXtalId(hit.volumeID()));
    
  m_calRelMap.insert(CalRelation(xtalIdx, &hit));

  return StatusCode::SUCCESS;
}

/// Inform that a new incident has occured
void CalSignalTool::handle ( const Incident& inc ) { 
  if ((inc.type() == "BeginEvent"))
    newEvent();

  return; 
}


/// apply electronic noise calcuation to entire cal
StatusCode CalSignalTool::calcNoise() {
  using namespace CalUtil;
  
  /// Loop through (installed) towers and crystals; 
  for (unsigned twrSeq = 0; twrSeq < m_twrList.size(); twrSeq++) {
    // get bay id of nth live tower
    const TwrNum twr(m_twrList[twrSeq]);
    for (LyrNum lyr; lyr.isValid(); lyr++)
      for (ColNum col; col.isValid(); col++) {
        
        // assemble current calXtalId
        const XtalIdx xtalIdx(twr.val(),
                              lyr.val(),
                              col.val());

        StatusCode sc = calcPoissonicNoiseXtal(xtalIdx);
        if (sc.isFailure())
          return sc;
        
        sc = calcElectronicNoiseXtal(xtalIdx);
        if (sc.isFailure())
          return sc;

        
      } // xtal loop
  } // twr loop
  
  return StatusCode::SUCCESS;
}

/// apply electronic noise calculation to single crystal
StatusCode CalSignalTool::calcElectronicNoiseXtal(const CalUtil::XtalIdx xtalIdx) {
  using namespace CalUtil;

  for (XtalDiode xtalDiode;
       xtalDiode.isValid();
       xtalDiode++) {

    const XtalRng xRng(xtalDiode.getFace(), 
                       xtalDiode.getDiode().getX8Rng());
    const RngIdx rngIdx(xtalIdx, xRng);
    
    // get pedestal sigma
    float pedSig;
    StatusCode sc = m_precalcCalib->getPedSigCIDAC(rngIdx, pedSig);
    if (sc.isFailure()) return sc;
    
    // use same rand for both ranges 
    // since the X1 & X8 noise
    const float rnd = CLHEP::RandGauss::shoot();
    
    m_calSignalMap[DiodeIdx(xtalIdx,xtalDiode)] += pedSig*rnd;
  } // for (XtalDiode)
  
  return StatusCode::SUCCESS;
}

/// apply poinssonic noise calculation to single cal crystal
StatusCode CalSignalTool::calcPoissonicNoiseXtal(const CalUtil::XtalIdx xtalIdx) {
  using namespace CalUtil;
  
  // retrieve mevPerDAC calibration object
  CalibData::CalMevPerDac const * const calibMPD = m_calCalibSvc->getMPD(xtalIdx);
  if (!calibMPD) return StatusCode::FAILURE;

  // store actual mevPerDAC values in usable array.
  CalArray<DiodeNum, float> mpd;
  mpd[CalUtil::LRG_DIODE] = calibMPD->getBig()->getVal();
  mpd[CalUtil::SM_DIODE]  = calibMPD->getSmall()->getVal();


  for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
    const DiodeNum diode = xDiode.getDiode();
    const DiodeIdx diodeIdx(xtalIdx, xDiode);
    
    // convert cidac in diode to MeV in xtal
    float meVXtal = m_calSignalMap[diodeIdx]*mpd[diode];
    
    // MeV in xtal -> electrons in diode
    float eDiode = meVXtal * m_ePerMeV[diode.val()];
    
    // apply poissonic fluctuation to # of electrons.
    const float noise = sqrt(eDiode)*CLHEP::RandGauss::shoot();
    
    // add noise
    eDiode += noise;
    
    // convert back to cidac in diode
    meVXtal = eDiode/m_ePerMeV[diode.val()];
    m_calSignalMap[diodeIdx] =
      meVXtal/mpd[diode];
    
  }
  return StatusCode::SUCCESS;
}


StatusCode CalSignalTool::getXtalSignalMap(const CalUtil::XtalIdx xtalIdx,
                                           XtalSignalMap &xtalSignalMap) {
  /// check that all internal data is up2date
  StatusCode sc = syncData();
  if (sc.isFailure())
    return 0;

  /// load up xtal array from full cal
  for (CalUtil::XtalDiode xtalDiode;
       xtalDiode.isValid();
       xtalDiode++)
    xtalSignalMap[xtalDiode] = m_calSignalMap[CalUtil::DiodeIdx(xtalIdx, xtalDiode)];
  
  return StatusCode::SUCCESS;
}
