// LOCAL
#include "PedMgr.h"

// GLAST
// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

/// retrieve pedestal vals for given xtal/face/rng
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

  CalibData::Ped *ped = getRangeBase(xtalId);

  //vals
  avr = ped->getAvr();
  sig = ped->getSig();
  cos = ped->getCosAngle();

  return StatusCode::SUCCESS;
}

StatusCode PedMgr::fillRangeBases() {
  m_rngBases.resize(RngIdx::N_VALS);

  for (RngIdx rngIdx; rngIdx.isValid(); rngIdx++) {
    CalXtalId xtalId = rngIdx.getCalXtalId();
    CalibData::RangeBase *rngBase = m_calibBase->getRange(xtalId);
    if (!rngBase) return StatusCode::FAILURE;

    if (!validateRangeBase(xtalId,rngBase)) return StatusCode::FAILURE;

    m_rngBases[rngIdx] = rngBase;
  }

  return StatusCode::SUCCESS;
}


StatusCode PedMgr::loadIdealVals() {
  MsgStream msglog(m_msgSvc, *m_logName);

  //-- SANITY CHECK --//
  if (m_idealCalib.pedVals.size() != RngNum::N_VALS) {
    msglog << MSG::ERROR << "wrong # of ideal pedestal vals." << endl;
    return StatusCode::FAILURE;;
  }
  if (m_idealCalib.pedCos.size() != RngNum::N_VALS) {
    msglog << MSG::ERROR << "wrong # of ideal ped cosine vals." << endl;
    return StatusCode::FAILURE;;
  }

  for (RngNum rng; rng.isValid(); rng++) {
    m_idealPeds[rng]   = m_idealCalib.pedVals[rng];
    m_idealPedSig[rng] = m_idealCalib.pedVals[rng] *
      m_idealCalib.pedSigPct;
    m_idealCos[rng] = m_idealCalib.pedCos[rng];
  }
  
  return StatusCode::SUCCESS;
}
