// LOCAL
#include "PedMgr.h"

// GLAST
// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

/// retrieve pedestal values for given xtal/face/rng
StatusCode PedMgr::getPed(const CalXtalId &xtalId,
                          float &avr,
                          float &sig,
                          float &cos) {

  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  CalibData::Ped *ped = getRangeBase(xtalId);
  if (!ped) return StatusCode::FAILURE;

  //values
  avr      = ped->getAvr();
  sig      = ped->getSig();
  cos      = ped->getCosAngle();

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
