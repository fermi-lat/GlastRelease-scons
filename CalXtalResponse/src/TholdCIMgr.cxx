// LOCAL
#include "TholdCIMgr.h"

// GLAST
// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

/// retrieve threshold calibration constants as measured w/ charnge injection
StatusCode TholdCIMgr::getTholds(const CalXtalId &xtalId,
                                 CalibData::ValSig &FLE,
                                 CalibData::ValSig &FHE,
                                 CalibData::ValSig &LAC
                                 ) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  if (m_idealMode) {
    FLE = m_idealFLE;
    FHE = m_idealFHE;
    LAC = m_idealLAC;
    return StatusCode::SUCCESS;
  }

  CalibData::CalTholdCI *tholdCI = getRangeBase(xtalId);
  if (!tholdCI) return StatusCode::FAILURE;

  //vals
  FLE = *(tholdCI->getFLE());
  FHE = *(tholdCI->getFHE());
  LAC = *(tholdCI->getLAC());

  return StatusCode::SUCCESS;
}

/// retrieve Upper Level Discriminator threshold as measured w/ charnge injection for given xtal/face/rng
StatusCode TholdCIMgr::getULD(const CalXtalId &xtalId,
                              CalibData::ValSig &ULDThold) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  if (m_idealMode) {
    ULDThold = m_idealULD[xtalId.getRange()];
    return StatusCode::SUCCESS;
  }

  CalibData::CalTholdCI *tholdCI = getRangeBase(xtalId);
  if (!tholdCI) return StatusCode::FAILURE;

  //vals
  ULDThold = *(tholdCI->getULD(xtalId.getRange()));

  return StatusCode::SUCCESS;
}

/// retrieve pedestal calibration constants as measured during charge injection threshold testing.
StatusCode TholdCIMgr::getPed(const CalXtalId &xtalId,
                              CalibData::ValSig &ped) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
  
  if (m_idealMode) {
    ped = m_idealPed[xtalId.getRange()];
    return StatusCode::SUCCESS;
  }

  CalibData::CalTholdCI *tholdCI = getRangeBase(xtalId);
  if (!tholdCI) return StatusCode::FAILURE;

  ped = *(tholdCI->getPed(xtalId.getRange()));

  return StatusCode::SUCCESS;
}

StatusCode TholdCIMgr::fillRangeBases() {
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

StatusCode TholdCIMgr::loadIdealVals() {
  MsgStream msglog(m_msgSvc, *m_logName); 
  
  //-- SANITY CHECKS --//
  if (m_idealCalib.ciULD.size() != RngNum::N_VALS) {
    msglog << MSG::ERROR << "Wrong # of ideal ULD vals." << endl;
    return StatusCode::FAILURE;;
  }
  if (m_idealCalib.ciPeds.size() != RngNum::N_VALS) {
    msglog << MSG::ERROR << "Wrong # of ideal ci Pedestal vals." << endl;
    return StatusCode::FAILURE;;
  }

  m_idealFLE.m_val = m_idealCalib.ciFLE;
  m_idealFLE.m_sig = m_idealCalib.ciFLE *
    m_idealCalib.ciSigPct;

  m_idealFHE.m_val = m_idealCalib.ciFHE;
  m_idealFHE.m_sig = m_idealCalib.ciFHE *
    m_idealCalib.ciSigPct;

  m_idealLAC.m_val = m_idealCalib.ciLAC;
  m_idealLAC.m_sig = m_idealCalib.ciLAC *
    m_idealCalib.ciSigPct;
  
  for (RngNum rng; rng.isValid(); rng++) {
    m_idealULD[rng].m_val = m_idealCalib.ciULD[rng];
    m_idealULD[rng].m_sig = m_idealCalib.ciULD[rng] *
      m_idealCalib.ciSigPct;

    m_idealPed[rng].m_val = m_idealCalib.ciPeds[rng];
    m_idealPed[rng].m_sig = m_idealCalib.ciPeds[rng] *
      m_idealCalib.ciSigPct;
  }

  return StatusCode::SUCCESS;
}
