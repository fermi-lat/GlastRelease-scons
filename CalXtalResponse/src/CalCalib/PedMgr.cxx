// LOCAL
#include "PedMgr.h"
#include "CalCalibSvc.h"

// GLAST
// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

/// get pedestal vals for given xtal/face/rng
StatusCode PedMgr::getPed(const CalXtalId &xtalId,
                          float &avr,
                          float &sig,
                          float &cos) {

  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
  
  if (m_idealMode) {
    RngNum rng = xtalId.getRange();
    // return default vals if we're in ideal (fake) mode
    avr = m_idealPeds[rng];
    sig = m_idealPedSig[rng];
    cos = m_idealCos[rng];
    return StatusCode::SUCCESS;
  } 

  CalibData::Ped *ped 
	  = (CalibData::Ped *)getRangeBase(xtalId);
  if (!ped) return StatusCode::FAILURE;

  //vals
  avr = ped->getAvr();
  sig = ped->getSig();
  cos = ped->getCosAngle();

  return StatusCode::SUCCESS;
}

StatusCode PedMgr::fillRangeBases() {
  m_rngBases.resize(RngIdx::N_VALS,0);

  for (RngIdx rngIdx; rngIdx.isValid(); rngIdx++) {
    CalXtalId xtalId = rngIdx.getCalXtalId();
    CalibData::RangeBase *rngBase = m_calibBase->getRange(xtalId);
    if (!rngBase) continue; // support partial LAT inst

    // support missing towers & missing crystals
    // keep moving if we're missing a particular calibration
    if (!validateRangeBase(rngBase)) continue;

    m_rngBases[rngIdx] = rngBase;
  }

  return StatusCode::SUCCESS;
}


StatusCode PedMgr::loadIdealVals() {

  //-- SANITY CHECK --//
  if (owner->m_idealCalib.pedVals.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name());
    msglog << MSG::ERROR << "wrong # of ideal pedestal vals." << endreq;
    return StatusCode::FAILURE;;
  }
  if (owner->m_idealCalib.pedCos.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name());
    msglog << MSG::ERROR << "wrong # of ideal ped cosine vals." << endreq;
    return StatusCode::FAILURE;;
  }
  if (owner->m_idealCalib.pedSigs.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name());
    msglog << MSG::ERROR << "wrong # of ideal ped sigma vals." << endreq;
    return StatusCode::FAILURE;;
  }

  for (RngNum rng; rng.isValid(); rng++) {
    m_idealPeds[rng]   = owner->m_idealCalib.pedVals[rng];
    m_idealPedSig[rng] = owner->m_idealCalib.pedSigs[rng];
    m_idealCos[rng]    = owner->m_idealCalib.pedCos[rng];
  }
  
  return StatusCode::SUCCESS;
}
