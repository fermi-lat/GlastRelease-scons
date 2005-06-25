// LOCAL
#include "TholdCIMgr.h"
#include "CalCalibSvc.h"

// GLAST
// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

/// get threshold calibration constants as measured w/ charnge injection
StatusCode TholdCIMgr::getTholds(const CalXtalId &xtalId,
                                 CalibData::ValSig &FLE,
                                 CalibData::ValSig &FHE,
                                 CalibData::ValSig &LAC
                                 ) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  if (m_idealMode) {
    FLE = m_idealFLE;
    FHE = m_idealFHE;
    LAC = m_idealLAC;
    return StatusCode::SUCCESS;
  }

  CalibData::CalTholdCI *tholdCI 
	  = (CalibData::CalTholdCI *)getRangeBase(xtalId);
  if (!tholdCI) return StatusCode::FAILURE;

  //vals
  FLE = *(tholdCI->getFLE());
  FHE = *(tholdCI->getFHE());
  LAC = *(tholdCI->getLAC());

  return StatusCode::SUCCESS;
}

/// get Upper Level Discriminator threshold as measured w/ charnge injection for given xtal/face/rng
StatusCode TholdCIMgr::getULD(const CalXtalId &xtalId,
                              CalibData::ValSig &ULDThold) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  if (m_idealMode) {
    ULDThold = m_idealULD[xtalId.getRange()];
    return StatusCode::SUCCESS;
  }

  CalibData::CalTholdCI *tholdCI 
	  = (CalibData::CalTholdCI *)getRangeBase(xtalId);
  if (!tholdCI) return StatusCode::FAILURE;

  //vals
  ULDThold = *(tholdCI->getULD(xtalId.getRange()));

  return StatusCode::SUCCESS;
}

/// get pedestal calibration constants as measured during charge injection threshold testing.
StatusCode TholdCIMgr::getPed(const CalXtalId &xtalId,
                              CalibData::ValSig &ped) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
  
  if (m_idealMode) {
    ped = m_idealPed[xtalId.getRange()];
    return StatusCode::SUCCESS;
  }

  CalibData::CalTholdCI *tholdCI 
	  = (CalibData::CalTholdCI *)getRangeBase(xtalId);
  if (!tholdCI) return StatusCode::FAILURE;

  ped = *(tholdCI->getPed(xtalId.getRange()));

  return StatusCode::SUCCESS;
}

StatusCode TholdCIMgr::fillRangeBases() {
  m_rngBases.resize(FaceIdx::N_VALS,0);

  for (FaceIdx faceIdx; faceIdx.isValid(); faceIdx++) {
    CalXtalId xtalId = faceIdx.getCalXtalId();
    CalibData::RangeBase *rngBase = m_calibBase->getRange(xtalId);
    if (!rngBase) continue; // support parial LAT inst

    // support missing towers & missing crystals
    // keep moving if we're missing a particular calibration
    if (!validateRangeBase(rngBase)) continue;

    m_rngBases[faceIdx] = rngBase;
  }

  return StatusCode::SUCCESS;
}

StatusCode TholdCIMgr::loadIdealVals() {
  
  //-- SANITY CHECKS --//
  if (owner->m_idealCalib.ciULD.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "Wrong # of ideal ULD vals." << endreq;
    return StatusCode::FAILURE;;
  }
  if (owner->m_idealCalib.ciPeds.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "Wrong # of ideal ci Pedestal vals." << endreq;
    return StatusCode::FAILURE;;
  }

  m_idealFLE.m_val = owner->m_idealCalib.ciFLE;
  m_idealFLE.m_sig = owner->m_idealCalib.ciFLE *
    owner->m_idealCalib.ciSigPct;

  m_idealFHE.m_val = owner->m_idealCalib.ciFHE;
  m_idealFHE.m_sig = owner->m_idealCalib.ciFHE *
    owner->m_idealCalib.ciSigPct;

  m_idealLAC.m_val = owner->m_idealCalib.ciLAC;
  m_idealLAC.m_sig = owner->m_idealCalib.ciLAC *
    owner->m_idealCalib.ciSigPct;
  
  for (RngNum rng; rng.isValid(); rng++) {
    m_idealULD[rng].m_val = owner->m_idealCalib.ciULD[rng];
    m_idealULD[rng].m_sig = owner->m_idealCalib.ciULD[rng] *
      owner->m_idealCalib.ciSigPct;

    m_idealPed[rng].m_val = owner->m_idealCalib.ciPeds[rng];
    m_idealPed[rng].m_sig = owner->m_idealCalib.ciPeds[rng] *
      owner->m_idealCalib.ciSigPct;
  }

  return StatusCode::SUCCESS;
}

bool TholdCIMgr::validateRangeBase(CalibData::RangeBase *rngBase) {
  CalibData::CalTholdCI *tholdCI = (CalibData::CalTholdCI*)(rngBase);

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
    MsgStream msglog(owner->msgSvc(), owner->name()); 
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

  if (peds->size() != (unsigned)RngNum::N_VALS ||
      ulds->size() != (unsigned) RngNum::N_VALS) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "can't get calib data for " 
           << m_calibPath;
	msglog << endreq;
    return false;
  }
  return true;
}
