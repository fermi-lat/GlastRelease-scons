// LOCAL
#include "TholdMuonMgr.h"

// GLAST
// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

bool TholdMuonMgr::validateRangeBase(const CalXtalId &xtalId, CalibData::RangeBase *rngBase) {
  // recast for specific calibration type
  CalibData::CalTholdMuon *tholdMuon = (CalibData::CalTholdMuon*)(rngBase);
  return true;
}

/// retrieve threshold calibration constants as measured w/ muon calibration
StatusCode TholdMuonMgr::getTholds(const CalXtalId &xtalId,
                                   CalibData::ValSig &FLE,
                                   CalibData::ValSig &FHE) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  CalibData::CalTholdMuon *tholdMuon = getRangeBase(xtalId);
  if (!tholdMuon) return StatusCode::FAILURE;

  //values
  FLE = *(tholdMuon->getFLE());
  FHE = *(tholdMuon->getFHE());

  return StatusCode::SUCCESS;
}

/// retrieve pedestal calibration constants as measured during muon calibration threshold testing.
StatusCode TholdMuonMgr::getPed(const CalXtalId &xtalId,
                                CalibData::ValSig &ped) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

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
