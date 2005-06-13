// LOCAL
#include "TholdMuonMgr.h"
#include "CalCalibSvc.h"

// GLAST
// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

/// get threshold calibration constants as measured w/ muon calibration
StatusCode TholdMuonMgr::getTholds(const CalXtalId &xtalId,
                                   CalibData::ValSig &FLE,
                                   CalibData::ValSig &FHE) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  if (m_idealMode) {
    FLE = m_idealFLE;
    FHE = m_idealFHE;
    return StatusCode::SUCCESS;
  }

  CalibData::CalTholdMuon *tholdMuon 
	  = (CalibData::CalTholdMuon *)getRangeBase(xtalId);
  if (!tholdMuon) return StatusCode::FAILURE;

  //vals
  FLE = *(tholdMuon->getFLE());
  FHE = *(tholdMuon->getFHE());

  return StatusCode::SUCCESS;
}

/** \brief get pedestal calibration constants as measured with muons
 */
StatusCode TholdMuonMgr::getPed(const CalXtalId &xtalId,
                                CalibData::ValSig &ped) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  if (m_idealMode) {
    ped = m_idealPed[xtalId.getRange()];
    return StatusCode::SUCCESS;
  }

  CalibData::CalTholdMuon *tholdMuon 
	  = (CalibData::CalTholdMuon *)getRangeBase(xtalId);
  if (!tholdMuon) return StatusCode::FAILURE;

  ped = *(tholdMuon->getPed(xtalId.getRange()));

  return StatusCode::SUCCESS;
}

StatusCode TholdMuonMgr::fillRangeBases() {
  m_rngBases.resize(FaceIdx::N_VALS,0);

  for (FaceIdx faceIdx; faceIdx.isValid(); faceIdx++) {
    CalXtalId xtalId = faceIdx.getCalXtalId();
    CalibData::RangeBase *rngBase = m_calibBase->getRange(xtalId);
    if (!rngBase) continue; // support partial LAT instruments

    // support missing towers & missing crystals
    // keep moving if we're missing a particular calibration
    if (!validateRangeBase(rngBase)) continue;

    m_rngBases[faceIdx] = rngBase;
  }

  return StatusCode::SUCCESS;
}

StatusCode TholdMuonMgr::loadIdealVals() {

  //-- SANITY CHECKS --//
  if (owner->m_idealCalib.muonPeds.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "Wrong # of ideal muon pedestal vals." << endreq;
    return StatusCode::FAILURE;;
  }

  m_idealFLE.m_val = owner->m_idealCalib.muonFLE;
  m_idealFLE.m_sig = owner->m_idealCalib.muonFLE *
    owner->m_idealCalib.muonSigPct;

  m_idealFHE.m_val = owner->m_idealCalib.muonFHE;
  m_idealFHE.m_sig = owner->m_idealCalib.muonFHE *
    owner->m_idealCalib.muonSigPct;
  
  for (RngNum rng; rng.isValid(); rng++) {
    m_idealPed[rng].m_val = owner->m_idealCalib.muonPeds[rng];
    m_idealPed[rng].m_sig = owner->m_idealCalib.muonPeds[rng] *
      owner->m_idealCalib.muonSigPct;
  }

  return StatusCode::SUCCESS;
}

bool TholdMuonMgr::validateRangeBase(CalibData::RangeBase *rngBase) {
  CalibData::CalTholdMuon *tholdMuon = (CalibData::CalTholdMuon*)(rngBase);

  if (tholdMuon->getFLE()) {
    // no error print out req'd b/c we're supporting LAT configs w/ empty bays
    // however, if tholdMuon->getFLE() is successful & following checks fail
    // then we have a problem b/c we have calib data which is only good for
    // partial xtal.
    return false;
  }
  if (tholdMuon->getFHE()) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "can't get calib data for " 
           << m_calibPath;
	msglog << endreq;
    return false;
  }
  
  const vector<ValSig> *peds = tholdMuon->getPeds();
  if (!peds) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "can't get calib data for " 
           << m_calibPath;
	msglog << endreq;
    return false;
  }
  if (peds->size() != RngNum::N_VALS) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "can't get calib data for " 
           << m_calibPath;
	msglog << endreq;
    return false;
  }

  return true;
}
