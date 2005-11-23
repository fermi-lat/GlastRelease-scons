// LOCAL
#include "TholdMuonMgr.h"
#include "CalCalibSvc.h"

// GLAST
// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

/// get threshold calibration constants as measured w/ muon calibration
StatusCode TholdMuonMgr::getTholds(CalXtalId xtalId,
                                   CalibData::ValSig &FLE,
                                   CalibData::ValSig &FHE) {
  // generic xtalId check
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
  
  // getTholds() specific xtalId check
  if (xtalId.validRange())
    throw invalid_argument("TholdMuon::getTholds() cannot accept range info in CalXtalId."
                           " Programmer error.");

  if (m_idealMode) {
    FLE = m_idealFLE;
    FHE = m_idealFHE;
    return StatusCode::SUCCESS;
  }

  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return sc;


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
StatusCode TholdMuonMgr::getPed(CalXtalId xtalId,
                                CalibData::ValSig &ped) {
  // generic xtalId check
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  // getPed() specific xtalId check
  if (!xtalId.validRange())
    throw invalid_argument("ThodMuon::getPed() requires range info in CalXtalId."
                          " Programmer error.");

  if (m_idealMode) {
    ped = m_idealPed[xtalId.getRange()];
    return StatusCode::SUCCESS;
  }

  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return sc;


  // need to create an xtalId object w/out range info...
  CalXtalId faceXtalId(xtalId.getTower(),
                       xtalId.getLayer(),
                       xtalId.getColumn(),
                       xtalId.getFace());
  CalibData::CalTholdMuon *tholdMuon 
    = (CalibData::CalTholdMuon *)getRangeBase(faceXtalId);
  if (!tholdMuon) return StatusCode::FAILURE;

  ped = *(tholdMuon->getPed(xtalId.getRange()));

  return StatusCode::SUCCESS;
}

StatusCode TholdMuonMgr::loadIdealVals() {

  //-- SANITY CHECKS --//
  if (owner->m_idealCalib.muonPeds.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
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


bool TholdMuonMgr::checkXtalId(CalXtalId xtalId) {
  // all queries require face info (leave range optional)
  if (!xtalId.validFace())
    throw invalid_argument("TholdMuon calib_type requires valid face info."
                           " Programmer error");
  return true;
}

bool TholdMuonMgr::validateRangeBase(CalibData::RangeBase *rangeBase) {
  CalibData::CalTholdMuon *tholdMuon = (CalibData::CalTholdMuon*)rangeBase;

  if (!tholdMuon->getFLE()) {
    // no error print out req'd b/c we're supporting LAT configs w/ empty bays
    // however, if tholdMuon->getFLE() is successful & following checks fail
    // then we have a problem b/c we have calib data which is only good for
    // partial xtal.
    return false;
  }
  if (!tholdMuon->getFHE()) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "can't get calib data for " 
           << m_calibPath;
    msglog << endreq;
    return false;
  }

  const vector<ValSig> *peds = tholdMuon->getPeds();
  if (!peds) {
    // no msg, b/c sometimes CalibSvc returns 'empty' TholdCI
    return false;
  }

  if (peds->size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "can't get calib data for " 
           << m_calibPath;
    msglog << endreq;
    return false;
  }
  return true;
}
