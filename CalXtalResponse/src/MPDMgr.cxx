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

/// retrieve MeVPerDac ratios for given xtal
StatusCode MPDMgr::getMPD(const CalXtalId &xtalId,
                          CalibData::ValSig &mpdLrg,
                          CalibData::ValSig &mpdSm) {

  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  if (m_idealMode) {
    mpdLrg = m_idealMPDLrg;
    mpdSm  = m_idealMPDSm;
    return StatusCode::SUCCESS;
  }

  // Retrieve generic pointer to rng specific data
  CalibData::CalMevPerDac *mpd = getRangeBase(xtalId);
  if (!mpd) return StatusCode::FAILURE;

  mpdLrg = *(mpd->getBig());
  mpdSm = *(mpd->getSmall());

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

StatusCode MPDMgr::loadIdealVals() {
  m_idealMPDLrg.m_val = m_idealCalib.mpdLrg;
  m_idealMPDLrg.m_sig = m_idealCalib.mpdLrg * m_idealCalib.mpdSigPct;

  m_idealMPDSm.m_val = m_idealCalib.mpdSm;
  m_idealMPDSm.m_sig = m_idealCalib.mpdSm * m_idealCalib.mpdSigPct;
  
  return StatusCode::SUCCESS;
}
