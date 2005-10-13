// LOCAL INCLUDES
#include "TholdCIMgr.h"
#include "CalCalibSvc.h"

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES

using namespace CalDefs;
using namespace idents;

/// get threshold calibration constants as measured w/ charnge injection
StatusCode TholdCIMgr::getTholds(CalXtalId xtalId,
                                 CalibData::ValSig &FLE,
                                 CalibData::ValSig &FHE,
                                 CalibData::ValSig &LAC
                                 ) {
  // generic xtalId check
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  // getTolds() specific xtalId check.
  if (xtalId.validRange())
    throw invalid_argument("TholdCI::getTholds() cannot accept range info in CalXtalId."
                           " Programmer error.");

  if (m_idealMode) {
    FLE = m_idealFLE;
    FHE = m_idealFHE;
    LAC = m_idealLAC;
    return StatusCode::SUCCESS;
  }

  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return sc;


  CalibData::CalTholdCI *tholdCI 
    = (CalibData::CalTholdCI *)getRangeBase(xtalId);
  if (!tholdCI) return StatusCode::FAILURE;

  //vals
  FLE = *(tholdCI->getFLE());
  FHE = *(tholdCI->getFHE());
  LAC = *(tholdCI->getLAC());

  return StatusCode::SUCCESS;
}

/// get Upper Level Discriminator threshold as measured w/ charnge injection for given xtal/face/rng
StatusCode TholdCIMgr::getULD(CalXtalId xtalId,
                              CalibData::ValSig &ULDThold) {
  // generic xtalId check
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  // getULD() specific xtalId check
  if (!xtalId.validRange())
    throw invalid_argument("ThodCI::getULD() requires range info in CalXtalId."
                          " Programmer error.");

  if (m_idealMode) {
    ULDThold = m_idealULD[xtalId.getRange()];
    return StatusCode::SUCCESS;
  }

  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return sc;

  // need to create an xtalId object w/out range info...
  CalXtalId faceXtalId(xtalId.getTower(),
                       xtalId.getLayer(),
                       xtalId.getColumn(),
                       xtalId.getFace());
  CalibData::CalTholdCI *tholdCI 
    = (CalibData::CalTholdCI *)getRangeBase(faceXtalId);
  if (!tholdCI) return StatusCode::FAILURE;

  //vals
  ULDThold = *(tholdCI->getULD(xtalId.getRange()));

  return StatusCode::SUCCESS;
}

/// get pedestal calibration constants as measured during charge injection threshold testing.
StatusCode TholdCIMgr::getPed(CalXtalId xtalId,
                              CalibData::ValSig &ped) {
  // generic xtalId check
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
  
  // getPed() specific xtalId check
  if (!xtalId.validRange())
    throw invalid_argument("ThodCI::getPed() requires range info in CalXtalId."
                          " Programmer error.");

  if (m_idealMode) {
    ped = m_idealPed[xtalId.getRange()];
    return StatusCode::SUCCESS;
  }

  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return sc;


  // need to create an xtalId object w/out range info...
  CalXtalId faceXtalId(xtalId.getTower(),
                       xtalId.getLayer(),
                       xtalId.getColumn(),
                       xtalId.getFace());
  CalibData::CalTholdCI *tholdCI 
    = (CalibData::CalTholdCI *)getRangeBase(faceXtalId);
  if (!tholdCI) return StatusCode::FAILURE;

  ped = *(tholdCI->getPed(xtalId.getRange()));

  return StatusCode::SUCCESS;
}

StatusCode TholdCIMgr::loadIdealVals() {
  
  //-- SANITY CHECKS --//
  if (owner->m_idealCalib.ciULD.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "Wrong # of ideal ULD vals." << endreq;
    return StatusCode::FAILURE;;
  }
  if (owner->m_idealCalib.ciPeds.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "Wrong # of ideal ci Pedestal vals." << endreq;
    return StatusCode::FAILURE;;
  }

  m_idealFLE.m_val = owner->m_idealCalib.ciFLE;
  m_idealFLE.m_sig = owner->m_idealCalib.ciFLE *
    owner->m_idealCalib.ciSigPct;

  m_idealFHE.m_val = owner->m_idealCalib.ciFHE;
  m_idealFHE.m_sig = owner->m_idealCalib.ciFHE *
    owner->m_idealCalib.ciSigPct;

  m_idealLAC.m_val = owner->m_idealCalib.ciLAC;
  m_idealLAC.m_sig = owner->m_idealCalib.ciLAC *
    owner->m_idealCalib.ciSigPct;
  
  for (RngNum rng; rng.isValid(); rng++) {
    m_idealULD[rng].m_val = owner->m_idealCalib.ciULD[rng];
    m_idealULD[rng].m_sig = owner->m_idealCalib.ciULD[rng] *
      owner->m_idealCalib.ciSigPct;

    m_idealPed[rng].m_val = owner->m_idealCalib.ciPeds[rng];
    m_idealPed[rng].m_sig = owner->m_idealCalib.ciPeds[rng] *
      owner->m_idealCalib.ciSigPct;
  }

  return StatusCode::SUCCESS;
}

bool TholdCIMgr::checkXtalId(CalXtalId xtalId) {
  // all queries require face info (leave range optional)
  if (!xtalId.validFace())
    throw invalid_argument("TholdCI calib_type requires valid face information in CalXtalId."
                           " Programmer error");
  return true;
}
