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

  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return sc;

  // Get generic pointer to rng specific data
  CalibData::CalMevPerDac *mpd = 
	  (CalibData::CalMevPerDac *)getRangeBase(xtalId);
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
