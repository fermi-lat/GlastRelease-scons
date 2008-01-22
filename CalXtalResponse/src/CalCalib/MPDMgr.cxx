// $Header$
/** @file
    @author Z.Fewtrell
*/
// LOCAL
#include "MPDMgr.h"

// GLAST

// EXTLIB

// STD

using namespace std;
using namespace CalUtil;
using namespace idents;
using namespace CalibData;


/// get MeVPerDac ratios for given xtal
const CalMevPerDac *MPDMgr::getMPD(XtalIdx xtalIdx) {
  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return NULL;

  // Get generic pointer to rng specific data
  return (CalMevPerDac*)m_rngBases[xtalIdx];
}

StatusCode MPDMgr::loadIdealVals() {
  const ValSig mpdLrg(m_ccsShared.m_idealCalib.mpdLrg,
                m_ccsShared.m_idealCalib.mpdLrg * 
                m_ccsShared.m_idealCalib.mpdSigPct);

  const ValSig mpdSm(m_ccsShared.m_idealCalib.mpdSm,
               m_ccsShared.m_idealCalib.mpdSm * 
               m_ccsShared.m_idealCalib.mpdSigPct);

  m_idealMPD.reset(new CalMevPerDac(&mpdLrg, &mpdSm));
  
  return StatusCode::SUCCESS;
}


StatusCode MPDMgr::genLocalStore() {
  for (XtalIdx xtalIdx; xtalIdx.isValid(); xtalIdx++) {
    if (!m_idealMode) {
      CalibData::CalMevPerDac const*const mpd = (CalMevPerDac*)getRangeBase(xtalIdx.getCalXtalId());
      if (!mpd) continue;
      if (!validateRangeBase(mpd)) continue;
      
      m_rngBases[xtalIdx] = mpd;
    } else m_rngBases[xtalIdx] = m_idealMPD.get();
  }
  

  return StatusCode::SUCCESS;
}

bool MPDMgr::validateRangeBase(const CalMevPerDac *mpd) {
  const ValSig *big = mpd->getBig();
  const ValSig *small = mpd->getSmall();
  if (!big || !small) return false;
  if (big->getVal() <= 0) return false;
  if (small->getVal() <=0) return false;

  return true;
}
