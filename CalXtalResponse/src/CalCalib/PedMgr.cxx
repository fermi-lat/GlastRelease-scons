// $Header$
/** @file
    @author Zach Fewtrell
 */
// LOCAL
#include "PedMgr.h"
#include "CalCalibSvc.h"

// GLAST
// EXTLIB
// STD

using namespace CalUtil;
using namespace idents;

/// get pedestal vals for given xtal/face/rng
StatusCode PedMgr::getPed(RngIdx rngIdx,
                          float &avr,
                          float &sig,
                          float &cos) {

  if (m_idealMode) {
    RngNum rng(rngIdx.getRng());
    // return default vals if we're in ideal (fake) mode
    avr = m_idealPeds[rng];
    sig = m_idealPedSig[rng];
    cos = m_idealCos[rng];
    return StatusCode::SUCCESS;
  } 

  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return sc;


  CalibData::Ped *ped 
    = (CalibData::Ped*)m_rngBases[rngIdx];
  if (!ped) return StatusCode::FAILURE;

  //vals
  avr = ped->getAvr();
  sig = ped->getSig();
  cos = ped->getCosAngle();

  return StatusCode::SUCCESS;
}

StatusCode PedMgr::loadIdealVals() {

  //-- SANITY CHECK --//
  if (owner->m_idealCalib.pedVals.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(owner->msgSvc(), owner->name());
    msglog << MSG::ERROR << "wrong # of ideal pedestal vals." << endreq;
    return StatusCode::FAILURE;;
  }
  if (owner->m_idealCalib.pedCos.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(owner->msgSvc(), owner->name());
    msglog << MSG::ERROR << "wrong # of ideal ped cosine vals." << endreq;
    return StatusCode::FAILURE;;
  }
  if (owner->m_idealCalib.pedSigs.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(owner->msgSvc(), owner->name());
    msglog << MSG::ERROR << "wrong # of ideal ped sigma vals." << endreq;
    return StatusCode::FAILURE;;
  }

  for (RngNum rng; rng.isValid(); rng++) {
    m_idealPeds[rng]   = owner->m_idealCalib.pedVals[rng.getInt()];
    m_idealPedSig[rng] = owner->m_idealCalib.pedSigs[rng.getInt()];
    m_idealCos[rng]    = owner->m_idealCalib.pedCos[rng.getInt()];
  }
  
  return StatusCode::SUCCESS;
}



StatusCode PedMgr::genLocalStore() {
  m_rngBases.resize(RngIdx::N_VALS,0);

  if (!m_idealMode) {
    for (RngIdx rngIdx; rngIdx.isValid(); rngIdx++) {
      CalibData::Ped *ped = (CalibData::Ped*)getRangeBase(rngIdx.getCalXtalId());
      if (!ped) continue;
      if (!validateRangeBase(ped)) continue;

      m_rngBases[rngIdx] = ped;
    }
  }

  return StatusCode::SUCCESS;
}

bool PedMgr::validateRangeBase(CalibData::Ped *ped) {
  if (!ped) return false;
  if (ped->getAvr() <= 0) return false;
  return true;
}
