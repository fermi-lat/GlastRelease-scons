// $Header$
/** @file
    @author Z.Fewtrell
 */
// LOCAL INCLUDES
#include "TholdCIMgr.h"

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES

using namespace CalUtil;
using namespace idents;

/// get threshold calibration constants as measured w/ charnge injection
const CalibData::CalTholdCI *TholdCIMgr::getTholdCI(const FaceIdx faceIdx) {
  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return NULL;

  return (CalibData::CalTholdCI*)m_rngBases[faceIdx];
}

StatusCode TholdCIMgr::loadIdealVals() {
  
  //-- SANITY CHECKS --//
  if (m_ccsShared.m_idealCalib.ciULD.size() != RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name()); 
    msglog << MSG::ERROR << "Wrong # of ideal ULD vals." << endreq;
    return StatusCode::FAILURE;;
  }
  if (m_ccsShared.m_idealCalib.ciPeds.size() != RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name()); 
    msglog << MSG::ERROR << "Wrong # of ideal ci Pedestal vals." << endreq;
    return StatusCode::FAILURE;;
  }

  CalibData::ValSig FLE, FHE, LAC;
  std::vector<CalibData::ValSig> ULD(4), Ped(4);

  FLE.m_val = m_ccsShared.m_idealCalib.ciFLE;
  FLE.m_sig = m_ccsShared.m_idealCalib.ciFLE *
    m_ccsShared.m_idealCalib.ciSigPct;

  FHE.m_val = m_ccsShared.m_idealCalib.ciFHE;
  FHE.m_sig = m_ccsShared.m_idealCalib.ciFHE *
    m_ccsShared.m_idealCalib.ciSigPct;

  LAC.m_val = m_ccsShared.m_idealCalib.ciLAC;
  LAC.m_sig = m_ccsShared.m_idealCalib.ciLAC *
    m_ccsShared.m_idealCalib.ciSigPct;
  
  for (RngNum rng; rng.isValid(); rng++) {
    ULD[rng.val()].m_val = m_ccsShared.m_idealCalib.ciULD[rng.val()];
    ULD[rng.val()].m_sig = m_ccsShared.m_idealCalib.ciULD[rng.val()] *
      m_ccsShared.m_idealCalib.ciSigPct;

    Ped[rng.val()].m_val = m_ccsShared.m_idealCalib.ciPeds[rng.val()];
    Ped[rng.val()].m_sig = m_ccsShared.m_idealCalib.ciPeds[rng.val()] *
      m_ccsShared.m_idealCalib.ciSigPct;
  }

  m_idealTholdCI.reset(new CalibData::CalTholdCI(&ULD, &FLE, &FHE, &LAC, &Ped));

  return StatusCode::SUCCESS;
}

bool TholdCIMgr::validateRangeBase(CalibData::CalTholdCI const*const tholdCI) {
  if (!tholdCI) return false;

  if (!tholdCI->getFLE()) {
    // no error print out req'd b/c we're supporting LAT configs w/ empty bays
    // however, if tholdCI->getFLE() is successful & following checks fail
    // then we have a problem b/c we have calib data which is only good for
    // partial xtal.
    return false;
  }
  if (!tholdCI->getFHE() ||
      !tholdCI->getLAC()) {
    // create MsgStream only when needed for performance
    MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name()); 
    msglog << MSG::ERROR << "can't get calib data for " 
           << m_calibPath;
    msglog << endreq;
    return false;
  }

  const vector<CalibData::ValSig> *peds = tholdCI->getPeds();
  const vector<CalibData::ValSig> *ulds = tholdCI->getULDs();
  if (!peds || !ulds) {
    // no msg, b/c sometimes CalibSvc returns 'empty' TholdCI
    return false;
  }

  if (peds->size() != RngNum::N_VALS ||
      ulds->size() != RngNum::N_VALS) {
    // create MsgStream only when needed for performance
    MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name()); 
    msglog << MSG::ERROR << "can't get calib data for " 
           << m_calibPath;
    msglog << endreq;
    return false;
  }
  return true;
}

StatusCode  TholdCIMgr::genLocalStore() {
  for (FaceIdx faceIdx; faceIdx.isValid(); faceIdx++) {
    if (!m_idealMode) {
      CalibData::CalTholdCI const *const  tholdCI = (CalibData::CalTholdCI*)getRangeBase(faceIdx.getCalXtalId());
      if (!tholdCI) continue;
      if (!validateRangeBase(tholdCI)) continue;
      
      m_rngBases[faceIdx] = tholdCI;
    } else m_rngBases[faceIdx] = m_idealTholdCI.get();
  } 

  return StatusCode::SUCCESS;
}

