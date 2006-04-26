// $Header$
/** @file
    @author Zach Fewtrell
 */
// LOCAL
#include "PedMgr.h"

// GLAST
// EXTLIB
// STD

using namespace CalUtil;
using namespace idents;

/// get pedestal vals for given xtal/face/rng
const Ped *PedMgr::getPed(RngIdx rngIdx) {
  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return NULL;

  return (Ped*)m_rngBases[rngIdx];
}

StatusCode PedMgr::loadIdealVals() {

  //-- SANITY CHECK --//
  if (m_ccsShared.m_idealCalib.pedVals.size() != RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name());
    msglog << MSG::ERROR << "wrong # of ideal pedestal vals." << endreq;
    return StatusCode::FAILURE;;
  }
  if (m_ccsShared.m_idealCalib.pedCos.size() != RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name());
    msglog << MSG::ERROR << "wrong # of ideal ped cosine vals." << endreq;
    return StatusCode::FAILURE;;
  }
  if (m_ccsShared.m_idealCalib.pedSigs.size() != RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name());
    msglog << MSG::ERROR << "wrong # of ideal ped sigma vals." << endreq;
    return StatusCode::FAILURE;;
  }

  // load ped data into ideal mode storage
  for (RngNum rng; rng.isValid(); rng++) {
    Ped newPed(m_ccsShared.m_idealCalib.pedVals[rng.val()],
               m_ccsShared.m_idealCalib.pedSigs[rng.val()]);
    m_idealPeds[rng] = newPed;
  }
  
  return StatusCode::SUCCESS;
}

StatusCode PedMgr::genLocalStore() {
  m_rngBases.resize(RngIdx::N_VALS,0);

  for (RngIdx rngIdx; rngIdx.isValid(); rngIdx++) {
    if (!m_idealMode) {
      Ped *ped = (Ped*)getRangeBase(rngIdx.getCalXtalId());
      if (!ped) continue;
      if (!validateRangeBase(ped)) continue;

      m_rngBases[rngIdx] = ped;
    } else m_rngBases[rngIdx] = &(m_idealPeds[rngIdx.getRng()]);
  } 

  return StatusCode::SUCCESS;
}

bool PedMgr::validateRangeBase(Ped *ped) {
  if (!ped) return false;
  if (ped->getAvr() <= 0) return false;
  return true;
}
