// LOCAL
#include "MPDMgr.h"

// GLAST
#include "CalibData/Cal/Xpos.h"

// EXTLIB

// STD
#include <algorithm>

using namespace std;
using namespace CalDefs;
using namespace idents;

bool MPDMgr::validateRangeBase(const CalXtalId &xtalId, CalibData::RangeBase *rngBase) {
  // recast for specific calibration type
  CalibData::CalMevPerDac* MeVPerDac = (CalibData::CalMevPerDac*)(rngBase);
  return true;
}

/// retrieve MeVPerDac ratios for given xtal
StatusCode MPDMgr::getMPD(const CalXtalId &xtalId,
                          CalibData::ValSig &lrg,
                          CalibData::ValSig &sm) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  // Retrieve generic pointer to rng specific data
  CalibData::CalMevPerDac *mpd = getRangeBase(xtalId);
  if (!mpd) return StatusCode::FAILURE;

  lrg = *(mpd->getBig());
  sm = *(mpd->getSmall());

  return StatusCode::SUCCESS;
}

StatusCode MPDMgr::fillRangeBases() {
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
