// LOCAL
#include "TholdCIMgr.h"

// GLAST
// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

/// retrieve threshold calibration constants as measured w/ charnge injection
StatusCode TholdCIMgr::getTholds(const CalXtalId &xtalId,
                                 CalibData::ValSig &FLE,
                                 CalibData::ValSig &FHE,
                                 CalibData::ValSig &LAC
                                 ) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  CalibData::CalTholdCI *tholdCI = getRangeBase(xtalId);
  if (!tholdCI) return StatusCode::FAILURE;

  //values
  FLE = *(tholdCI->getFLE());
  FHE = *(tholdCI->getFHE());
  LAC = *(tholdCI->getLAC());

  return StatusCode::SUCCESS;
}

/// retrieve Upper Level Discriminator threshold as measured w/ charnge injection for given xtal/face/rng
StatusCode TholdCIMgr::getULD(const CalXtalId &xtalId,
                              CalibData::ValSig &ULDThold) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  CalibData::CalTholdCI *tholdCI = getRangeBase(xtalId);
  if (!tholdCI) return StatusCode::FAILURE;

  //values
  ULDThold = *(tholdCI->getULD(xtalId.getRange()));

  return StatusCode::SUCCESS;
}

/// retrieve pedestal calibration constants as measured during charge injection threshold testing.
StatusCode TholdCIMgr::getPed(const CalXtalId &xtalId,
                              CalibData::ValSig &ped) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  CalibData::CalTholdCI *tholdCI = getRangeBase(xtalId);
  if (!tholdCI) return StatusCode::FAILURE;

  ped = *(tholdCI->getPed(xtalId.getRange()));

  return StatusCode::SUCCESS;
}

StatusCode TholdCIMgr::fillRangeBases() {
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
