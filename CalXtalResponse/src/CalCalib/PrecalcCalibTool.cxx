// $Header$
/** @file
    @author Z.Fewtrell
*/


// LOCAL INCLUDES
#include "PrecalcCalibTool.h"
#include "CalCalibSvc.h"

// GLAST INCLUDES

// EXTLIB INCLUDES
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IIncidentSvc.h"

// STD INCLUDES
#include <algorithm>

using namespace CalUtil;
using namespace std;

//static ToolFactory<PrecalcCalibTool> s_factory;
//const IToolFactory& PrecalcCalibToolFactory = s_factory;
DECLARE_TOOL_FACTORY(PrecalcCalibTool);

/// used to represent invalid values in internal arrays.
const float BAD_FLOAT = -999999.999999F;

PrecalcCalibTool::PrecalcCalibTool( const string& type,
                                    const string& name,
                                    const IInterface* parent)
  : AlgTool(type,name,parent),
    m_serNoPed(-1),
    m_serNoINL(-1),
    m_serNoTholdCI(-1),
    m_initialized(false),
    m_isValid(false),
    m_calCalibSvc(0)

{
  declareInterface<IPrecalcCalibTool>(this);
  
  declareProperty("CalCalibSvc", m_calCalibSvcName = "CalCalibSvc");

  // initialize internal arrays
  using namespace std;
  fill(m_trigCIDAC.begin(),
       m_trigCIDAC.end(),
       BAD_FLOAT);
  fill(m_trigADC.begin(),
        m_trigADC.end(),
        BAD_FLOAT);
  fill(m_trigMeV.begin(),
       m_trigMeV.end(),
       BAD_FLOAT);
  fill(m_trigRng.begin(),
       m_trigRng.end(),
       LEX8);
  fill(m_lacCIDAC.begin(),
       m_lacCIDAC.end(),
       BAD_FLOAT);
  fill(m_pedSigCIDAC.begin(),
       m_pedSigCIDAC.end(),
       BAD_FLOAT);
}

/// intialize / retrieve all needed Gaudi based objects
StatusCode PrecalcCalibTool::initialize() {
  MsgStream msglog(msgSvc(), name());   
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  //-- jobOptions --//
  if ((sc = setProperties()).isFailure()) {
    msglog << MSG::ERROR << "Failed to set properties" << endreq;
    return sc;
  }

  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get " << m_calCalibSvcName << endreq;
    return sc;
  }

  // Get ready to listen for BeginEvent
  IIncidentSvc* incSvc;
  sc = service("IncidentSvc", incSvc, true);
  if (sc.isSuccess() ) {
    int priority = CalCalibSvc::INCIDENT_PRIORITY+1;  // this should be lower priority (higher #?) than CalCalibSvc
    incSvc->addListener(this, "BeginEvent", priority);
  } else {
    msglog << MSG::ERROR << "can't find IncidentSvc" << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}



/// check if any calibrations have changed out from under us
StatusCode PrecalcCalibTool::updateCalib() {
  StatusCode sc;

  // quickly return if we're already checked for this event.
  if (m_isValid) return StatusCode::SUCCESS;

  bool needUpdate = false;

  int tmpSerNo = m_calCalibSvc->getSerNoPed();
  if (m_serNoPed != tmpSerNo) needUpdate = true;
  m_serNoPed = tmpSerNo;

  tmpSerNo = m_calCalibSvc->getSerNoINL();
  if (m_serNoINL != tmpSerNo) needUpdate = true;
  m_serNoINL = tmpSerNo;

  tmpSerNo = m_calCalibSvc->getSerNoTholdCI();
  if (m_serNoTholdCI != tmpSerNo) needUpdate = true;
  m_serNoTholdCI = tmpSerNo;

  if (needUpdate || !m_initialized) {
    sc = genLocalStore();
    if (sc.isFailure()) return sc;

    m_initialized = true;
  }
  
  m_isValid = true;
  return StatusCode::SUCCESS;
}

StatusCode PrecalcCalibTool::genLocalStore() {
  StatusCode sc;

  for (CalUtil::FaceIdx faceIdx; faceIdx.isValid(); faceIdx++) {

    //-- TholdCI --//
    CalibData::CalTholdCI const*const tholdCI = m_calCalibSvc->getTholdCI(faceIdx);
    // allow for missing tholdci data, 1st time only (all towers may not be populated)
    if (!tholdCI) continue;

    //-- LAC--//
    sc = m_calCalibSvc->evalCIDAC(RngIdx(faceIdx, CalUtil::LEX8), 
                                  tholdCI->getLAC()->getVal(), 
                                  m_lacCIDAC[faceIdx]);
    // allow for missing intlin data, 1st time only (all towers may not be populated)
    if (sc.isFailure()) continue;

    //-- FLE --//
    // thresholds are in X8 units, but may need conversion to X1 units if they
    // are past hardware range.
    DiodeIdx diodeIdx(faceIdx, LRG_DIODE);
    m_trigADC[diodeIdx] = tholdCI->getFLE()->getVal();
    m_trigRng[diodeIdx] = (m_trigADC[diodeIdx] > tholdCI->getULD(LEX8.val())->getVal()) ?
      LEX1 : LEX8;
    RngIdx rngIdx(faceIdx, m_trigRng[diodeIdx]);
    
    // convert to LEX1 range if needed
    if (m_trigRng[diodeIdx] == LEX1) {
      sc = lex8_to_lex1(faceIdx, 
                        m_trigADC[diodeIdx], 
                        m_trigADC[diodeIdx]);
      // allow for missing lex8 to lex1 data, 1st time only (all towers may not be populated)
      if (sc.isFailure()) continue;
    }
    sc = m_calCalibSvc->evalCIDAC(rngIdx, 
                                  m_trigADC[diodeIdx], 
                                  m_trigCIDAC[diodeIdx]);
    // intlin should work if we got this far
    if (sc.isFailure()) return sc;

    // calculate face signal mev for trigger
    sc = m_calCalibSvc->evalFaceSignal(rngIdx,
                                       m_trigADC[diodeIdx],
                                       m_trigMeV[diodeIdx]);
    // allow for missing faceSignal data, 1st time only (all towers may not be populated)
    if (sc.isFailure()) continue;
    


    //-- FHE --//
    // thresholds are in X8 units, but may need conversion to X1 units if they
    // are past hardware range.
    diodeIdx = DiodeIdx(faceIdx, SM_DIODE);
    m_trigADC[diodeIdx] = tholdCI->getFHE()->getVal();
    m_trigRng[diodeIdx] = (m_trigADC[diodeIdx] > tholdCI->getULD(HEX8.val())->getVal()) ?
      HEX1 : HEX8;
    rngIdx = RngIdx(faceIdx, m_trigRng[diodeIdx]);
    
    // convert to HEX1 range if needed
    if (m_trigRng[diodeIdx] == HEX1) {
      sc = hex8_to_hex1(faceIdx,                         
                        m_trigADC[diodeIdx], 
                        m_trigADC[diodeIdx]);
      // allow for missing hex8_to_hex1 data, 1st time only (all towers may not be populated)
      if (sc.isFailure()) continue;
    }
    sc = m_calCalibSvc->evalCIDAC(rngIdx, 
                                  m_trigADC[diodeIdx], 
                                  m_trigCIDAC[diodeIdx]);
    // this should work by now if we got this far.
    if (sc.isFailure()) return sc;

    // calculate face signal mev for trigger
    sc = m_calCalibSvc->evalFaceSignal(rngIdx,
                                       m_trigADC[diodeIdx],
                                       m_trigMeV[diodeIdx]);
    // this should work by now if we got this far.
    if (sc.isFailure()) return sc;
    

    //-- (per range) Ped, ULDCI & PEDCI --//
    for (RngNum rng; rng.isValid(); rng++) {
      const RngIdx rngIdx(faceIdx, rng);

      //  CalibData::Ped const*const ped = m_calCalibSvc->getPed(rngIdx);
      float pedsig;
      sc = m_calCalibSvc->getPedSig(rngIdx,pedsig);
      if(sc.isFailure())return sc;
      
      float sigDAC;

      // this should work by now if we got this far.
      sc = m_calCalibSvc->evalCIDAC(rngIdx, pedsig, sigDAC);
      if (sc.isFailure()) return sc;

      // DAC scale is ped subtracted, so we'll make pedstal @ 0
      // ignore cos.
      m_pedSigCIDAC[rngIdx] = sigDAC;
    }  
  }
  
  return StatusCode::SUCCESS;
}


/// return pedestal sigma converted to CIDAC scale
StatusCode PrecalcCalibTool::getPedSigCIDAC(RngIdx rngIdx, float &pedSigCIDAC) {
  StatusCode sc;

  sc = updateCalib();
  if (sc.isFailure()) return sc;

  pedSigCIDAC = m_pedSigCIDAC[rngIdx];
  if (pedSigCIDAC == BAD_FLOAT) return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}

/// return trigger threshold in CIDAC scale
StatusCode PrecalcCalibTool::getTrigCIDAC(DiodeIdx diodeIdx, float &trigCIDAC) {
  StatusCode sc;

  sc = updateCalib();
  if (sc.isFailure()) return sc;
  
  trigCIDAC = m_trigCIDAC[diodeIdx];
  if (trigCIDAC == BAD_FLOAT) return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}

/// return trigger threshold in CIDAC scale
StatusCode PrecalcCalibTool::getTrigMeV(DiodeIdx diodeIdx, float &mev) {
  StatusCode sc;

  sc = updateCalib();
  if (sc.isFailure()) return sc;
  
  mev = m_trigMeV[diodeIdx];
  if (mev == BAD_FLOAT) return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}

/// return trigger threshold in CIDAC scale
StatusCode PrecalcCalibTool::getTrigRngADC(DiodeIdx diodeIdx, RngNum &rng, float &adc) {
  StatusCode sc;

  sc = updateCalib();
  if (sc.isFailure()) return sc;
  
  adc = m_trigADC[diodeIdx];
  if (adc == BAD_FLOAT) return StatusCode::FAILURE;
  rng = m_trigRng[diodeIdx];

  return StatusCode::SUCCESS;
}

/// return lac threshold in CIDAC scale
StatusCode PrecalcCalibTool::getLacCIDAC(FaceIdx faceIdx, float &lacCIDAC) {
  StatusCode sc;

  sc = updateCalib();
  if (sc.isFailure()) return sc;
  
  lacCIDAC = m_lacCIDAC[faceIdx];
  if (lacCIDAC == BAD_FLOAT) return StatusCode::FAILURE;
  
  return StatusCode::SUCCESS;
}

StatusCode PrecalcCalibTool::lex8_to_lex1(const FaceIdx faceIdx, 
                                          const float l8adc, 
                                          float &l1adc) {

  float tmpCIDAC;
  StatusCode sc;

  //-- TholdCI --//
  CalibData::CalTholdCI const*const tholdCI = m_calCalibSvc->getTholdCI(faceIdx);
  if (!tholdCI) return StatusCode::SUCCESS;
  
  // use this point to generate LEX8/LEX1 ratio
  const float x8tmp = tholdCI->getULD(LEX8.val())->getVal();

  // 1st convert to cidac
  sc = m_calCalibSvc->evalCIDAC(RngIdx(faceIdx, LEX8),
                                x8tmp, tmpCIDAC);
  if (sc.isFailure()) return sc;
      
  // 2nd convert to LEX1 adc
  float x1tmp;
  sc = m_calCalibSvc->evalADC(RngIdx(faceIdx, LEX1),
                              tmpCIDAC, x1tmp);
  if (sc.isFailure()) return sc;
      
  const float rat = x1tmp/x8tmp;

  l1adc = l8adc*rat;

  return StatusCode::SUCCESS;

}

StatusCode PrecalcCalibTool::hex8_to_hex1(const FaceIdx faceIdx, 
                                          const float h8adc, 
                                          float &h1adc) {

  float tmpCIDAC;
  StatusCode sc;

  //-- TholdCI --//
  CalibData::CalTholdCI const*const tholdCI = m_calCalibSvc->getTholdCI(faceIdx);
  if (!tholdCI) return StatusCode::SUCCESS;

  // use this point to generate HEX8/HEX1 ratio
  const float x8tmp = tholdCI->getULD(HEX8.val())->getVal();

  // 1st convert to cidac
  sc = m_calCalibSvc->evalCIDAC(RngIdx(faceIdx, HEX8),
                                x8tmp, tmpCIDAC);
  if (sc.isFailure()) return sc;
      
  // 2nd convert to LEX1 adc
  float x1tmp;
  sc = m_calCalibSvc->evalADC(RngIdx(faceIdx, HEX1),
                              tmpCIDAC, x1tmp);
  if (sc.isFailure()) return sc;
      
  const float rat = x1tmp/x8tmp;

  h1adc = h8adc*rat;

  return StatusCode::SUCCESS;

}


