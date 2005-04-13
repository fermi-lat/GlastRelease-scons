// LOCAL
#include "MPDMgr.h"
#include "CalCalibSvc.h"

// GLAST
#include "CalibData/Cal/Xpos.h"

// EXTLIB

// STD
#include <algorithm>

using namespace std;
using namespace CalDefs;
using namespace idents;

/// get MeVPerDac ratios for given xtal
StatusCode MPDMgr::getMPD(const CalXtalId &xtalId,
                          CalibData::ValSig &mpdLrg,
                          CalibData::ValSig &mpdSm) {

  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  if (m_idealMode) {
    mpdLrg = m_idealMPDLrg;
    mpdSm  = m_idealMPDSm;
    return StatusCode::SUCCESS;
  }

  // Get generic pointer to rng specific data
  CalibData::CalMevPerDac *mpd = getRangeBase(xtalId);
  if (!mpd) return StatusCode::FAILURE;

  mpdLrg = *(mpd->getBig());
  mpdSm = *(mpd->getSmall());

  return StatusCode::SUCCESS;
}

StatusCode MPDMgr::fillRangeBases() {
  m_rngBases.resize(XtalIdx::N_VALS,0);

  for (XtalIdx xtalIdx; xtalIdx.isValid(); xtalIdx++) {
    CalXtalId xtalId = xtalIdx.getCalXtalId();
    CalibData::RangeBase *rngBase = m_calibBase->getRange(xtalId);
    if (!rngBase) continue; // support partial LAT inst

    // support missing towers & missing crystals
    // keep moving if we're missing a particular calibration
    if (!validateRangeBase(xtalId,rngBase)) continue;

    m_rngBases[xtalIdx] = rngBase;
  }

  return StatusCode::SUCCESS;
}

StatusCode MPDMgr::loadIdealVals() {
  m_idealMPDLrg.m_val = owner->m_idealCalib.mpdLrg;
  m_idealMPDLrg.m_sig = owner->m_idealCalib.mpdLrg * 
     owner->m_idealCalib.mpdSigPct;

  m_idealMPDSm.m_val = owner->m_idealCalib.mpdSm;
  m_idealMPDSm.m_sig = owner->m_idealCalib.mpdSm * 
     owner->m_idealCalib.mpdSigPct;
  
  return StatusCode::SUCCESS;
}
