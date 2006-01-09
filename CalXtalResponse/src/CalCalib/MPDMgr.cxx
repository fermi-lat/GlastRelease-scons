// $Header$
/** @file
    @author Zach Fewtrell
 */
// LOCAL
#include "MPDMgr.h"
#include "CalCalibSvc.h"

// GLAST
#include "CalibData/Cal/Xpos.h"

// EXTLIB

// STD

using namespace std;
using namespace CalUtil;
using namespace idents;

/// get MeVPerDac ratios for given xtal
StatusCode MPDMgr::getMPD(XtalIdx xtalIdx,
                          CalibData::ValSig &mpdLrg,
                          CalibData::ValSig &mpdSm) {

  if (m_idealMode) {
    mpdLrg = m_idealMPDLrg;
    mpdSm  = m_idealMPDSm;
    return StatusCode::SUCCESS;
  }

  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return sc;

  // Get generic pointer to rng specific data
  CalibData::CalMevPerDac *mpd = 
    (CalibData::CalMevPerDac*)m_rngBases[xtalIdx];
  if (!mpd) return StatusCode::FAILURE;

  mpdLrg = *(mpd->getBig());
  mpdSm = *(mpd->getSmall());

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


StatusCode MPDMgr::genLocalStore() {
  m_rngBases.resize(XtalIdx::N_VALS,0);

  if (!m_idealMode) {
    for (XtalIdx xtalIdx; xtalIdx.isValid(); xtalIdx++) {
      CalibData::CalMevPerDac *mpd = (CalibData::CalMevPerDac*)getRangeBase(xtalIdx.getCalXtalId());
      if (!mpd) continue;
      if (!validateRangeBase(mpd)) continue;
      
      m_rngBases[xtalIdx] = mpd;
    }
  }
  return StatusCode::SUCCESS;
}

bool MPDMgr::validateRangeBase(CalibData::CalMevPerDac *mpd) {
  const ValSig *big = mpd->getBig();
  const ValSig *small = mpd->getSmall();
  if (!big || !small) return false;
  if (big->getVal() <= 0) return false;
  if (small->getVal() <=0) return false;

  return true;
}
