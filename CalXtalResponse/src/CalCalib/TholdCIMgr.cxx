// $Header$
/** @file
    @author Zach Fewtrell
 */
// LOCAL INCLUDES
#include "TholdCIMgr.h"

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES

using namespace CalUtil;
using namespace idents;

/// get threshold calibration constants as measured w/ charnge injection
const CalTholdCI *TholdCIMgr::getTholdCI(FaceIdx faceIdx) {
  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return NULL;

  return (CalTholdCI*)m_rngBases[faceIdx];
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

  ValSig FLE, FHE, LAC;
  vector<ValSig> ULD(4), Ped(4);

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

  m_idealTholdCI.reset(new CalTholdCI(&ULD, &FLE, &FHE, &LAC, &Ped));

  return StatusCode::SUCCESS;
}

bool TholdCIMgr::validateRangeBase(CalTholdCI *tholdCI) {
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

  const vector<ValSig> *peds = tholdCI->getPeds();
  const vector<ValSig> *ulds = tholdCI->getULDs();
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
  m_rngBases.resize(FaceIdx::N_VALS,0);

  for (FaceIdx faceIdx; faceIdx.isValid(); faceIdx++) {
    if (!m_idealMode) {
      CalTholdCI *tholdCI = (CalTholdCI*)getRangeBase(faceIdx.getCalXtalId());
      if (!tholdCI) continue;
      if (!validateRangeBase(tholdCI)) continue;
      
      m_rngBases[faceIdx] = tholdCI;
    } else m_rngBases[faceIdx] = m_idealTholdCI.get();
  } 

  return StatusCode::SUCCESS;
}

