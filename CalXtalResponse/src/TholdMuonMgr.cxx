// LOCAL
#include "TholdMuonMgr.h"

// GLAST
// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

/// retrieve threshold calibration constants as measured w/ muon calibration
StatusCode TholdMuonMgr::getTholds(const CalXtalId &xtalId,
                                   CalibData::ValSig &FLE,
                                   CalibData::ValSig &FHE) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  if (m_idealMode) {
    FLE = m_idealFLE;
    FHE = m_idealFHE;
    return StatusCode::SUCCESS;
  }

  CalibData::CalTholdMuon *tholdMuon = getRangeBase(xtalId);
  if (!tholdMuon) return StatusCode::FAILURE;

  //vals
  FLE = *(tholdMuon->getFLE());
  FHE = *(tholdMuon->getFHE());

  return StatusCode::SUCCESS;
}

/// retrieve pedestal calibration constants as measured during muon calibration threshold testing.
StatusCode TholdMuonMgr::getPed(const CalXtalId &xtalId,
                                CalibData::ValSig &ped) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  if (m_idealMode) {
    ped = m_idealPed[xtalId.getRange()];
    return StatusCode::SUCCESS;
  }

  CalibData::CalTholdMuon *tholdMuon = getRangeBase(xtalId);
  if (!tholdMuon) return StatusCode::FAILURE;

  ped = *(tholdMuon->getPed(xtalId.getRange()));

  return StatusCode::SUCCESS;
}

StatusCode TholdMuonMgr::fillRangeBases() {
  m_rngBases.resize(XtalIdx::N_VALS);

  for (XtalIdx xtalIdx; xtalIdx.isValid(); xtalIdx++) {
    CalXtalId xtalId = xtalIdx.getCalXtalId();
    CalibData::RangeBase *rngBase = m_calibBase->getRange(xtalId);
    if (!rngBase) return StatusCode::FAILURE;

    if (!validateRangeBase(xtalId,rngBase)) return StatusCode::FAILURE;

    m_rngBases[xtalIdx] = rngBase;
  }

  return StatusCode::SUCCESS;
}

StatusCode TholdMuonMgr::loadIdealVals() {
  MsgStream msglog(m_msgSvc, *m_logName); 

  //-- SANITY CHECKS --//
  if (m_idealCalib.muonPeds.size() != RngNum::N_VALS) {
    msglog << MSG::ERROR << "Wrong # of ideal muon pedestal vals." << endl;
    return StatusCode::FAILURE;;
  }

  m_idealFLE.m_val = m_idealCalib.muonFLE;
  m_idealFLE.m_sig = m_idealCalib.muonFLE *
    m_idealCalib.muonSigPct;

  m_idealFHE.m_val = m_idealCalib.muonFHE;
  m_idealFHE.m_sig = m_idealCalib.muonFHE *
    m_idealCalib.muonSigPct;
  
  for (RngNum rng; rng.isValid(); rng++) {
    m_idealPed[rng].m_val = m_idealCalib.muonPeds[rng];
    m_idealPed[rng].m_sig = m_idealCalib.muonPeds[rng] *
      m_idealCalib.muonSigPct;
  }

  return StatusCode::SUCCESS;
}
