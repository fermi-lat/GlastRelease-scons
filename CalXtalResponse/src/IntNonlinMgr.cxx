// LOCAL
#include "IntNonlinMgr.h"
#include "CalCalibSvc.h"

// GLAST
#include "CalibData/DacCol.h"

// EXTLIB

// STD
#include <algorithm>

using namespace CalDefs;
using namespace idents;

IntNonlinMgr::IntNonlinMgr() : 
  CalibItemMgr(CalibData::CAL_IntNonlin, 
               N_SPLINE_TYPES),
  m_idealADCs(RngNum::N_VALS), // one spline per range
  m_idealDACs(RngNum::N_VALS)  // one spline per range
{

  // set size of spline lists (1 per range)
  for (unsigned i = 0; i < m_splineLists.size(); i++) {
    m_splineLists[i].resize(RngIdx::N_VALS, 0);
    m_splineXMin[i].resize (RngIdx::N_VALS, 0);
    m_splineXMax[i].resize (RngIdx::N_VALS, 0);
  }
}

bool IntNonlinMgr::validateRangeBase(CalibData::RangeBase *rngBase) {
  // recast for specific calibration type
  CalibData::IntNonlin* intNonlin = (CalibData::IntNonlin*)(rngBase);

  //get vector of vals
  const vector<float> *intNonlinVec = intNonlin->getValues();
  if (!intNonlinVec)    
    return false;

  return true;
}

/// get integral non-linearity vals for given xtal/face/rng
StatusCode IntNonlinMgr::getIntNonlin(const CalXtalId &xtalId,
                                      const vector< float > *&adcs,
                                      const vector< unsigned > *&dacs,
                                      float &error) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  if (m_idealMode) {
    RngNum rng = xtalId.getRange();
    adcs = &m_idealADCs[rng];
    dacs = &m_idealDACs[rng];
    error = m_idealErr;
    return StatusCode::SUCCESS;
  }

  
  CalibData::IntNonlin *intNonlin 
	  = (CalibData::IntNonlin *)getRangeBase(xtalId);
  if (!intNonlin) return StatusCode::FAILURE;

  //get vector of vals
  const vector<float> *intNonlinVec = intNonlin->getValues();

  //get collection of associated DAC vals
  CalibData::DacCol *intNonlinDacCol = m_calibBase->getDacCol(xtalId.getRange());
  const vector<unsigned> *intNonlinDacVec;
  intNonlinDacVec = intNonlinDacCol->getDacs();

  // check array lens (can't check during validateRangeBase() since DAC is not
  // part of rangebase
  if (intNonlinVec->size() > intNonlinDacVec->size()) {
    return StatusCode::FAILURE;
  }

  //whew, we made it this far... let's assign our outputs & leave!
  adcs = intNonlinVec;
  dacs = intNonlinDacVec;
  error = intNonlin->getError();

  return StatusCode::SUCCESS;
}

StatusCode IntNonlinMgr::genSplines() {
  StatusCode sc;
  const vector<float> *adc;
  const vector<unsigned> *dac;
  vector<double> dblDac;
  vector<double> dblAdc;

  for (RngIdx rngIdx; rngIdx.isValid(); rngIdx++) {
    float error;
    
    // support missing towers & missing crystals
    // keep moving if we're missing a particular calibration
    sc = getIntNonlin(rngIdx.getCalXtalId(), adc, dac, error);
    if (sc.isFailure()) continue;

    int n = min(adc->size(),dac->size());

    dblDac.resize(n);
    dblAdc.resize(n);

    // create double vector for input to genSpline()
    copy(adc->begin(), adc->begin() + n, dblAdc.begin());
    copy(dac->begin(), dac->begin() + n, dblDac.begin());

	CalXtalId xtalId = rngIdx.getCalXtalId();

    // put rng id string into spline name
    ostringstream rngStr;
    rngStr << '[' << xtalId
           << ']';

    genSpline(INL_SPLINE, xtalId, "INL"    + rngStr.str(),    
              dblAdc, dblDac);
    genSpline(INV_INL_SPLINE, xtalId, "invINL" + rngStr.str(), 
              dblDac, dblAdc);
  }

  return StatusCode::SUCCESS;
}

StatusCode IntNonlinMgr::fillRangeBases() {
  m_rngBases.resize(RngIdx::N_VALS,0);

  for (RngIdx rngIdx; rngIdx.isValid(); rngIdx++) {
    CalXtalId xtalId = rngIdx.getCalXtalId();
    CalibData::RangeBase *rngBase = m_calibBase->getRange(xtalId);
    if (!rngBase) continue; // support partial LAT inst

    // support missing towers & missing crystals
    // keep moving if we're missing a particular calibration
    if (!validateRangeBase(rngBase)) continue;

    const vector<float> *adc = ((CalibData::IntNonlin*)rngBase)->getValues();

    m_rngBases[rngIdx] = rngBase;
  }

  return StatusCode::SUCCESS;
}

StatusCode IntNonlinMgr::loadIdealVals() {
  const int maxADC = 4095;

  //-- SANITY CHECKS --//
  if (owner->m_idealCalib.ciULD.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "Bad # of ULD vals in ideal CalCalib xml file" 
           << endreq;
    return StatusCode::FAILURE;
  }
  if (owner->m_idealCalib.inlADCPerDAC.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "Bad # of ADCPerDAC vals in ideal CalCalib xml file" 
           << endreq;
    return StatusCode::FAILURE;
  }
  
  for (RngNum rng; rng.isValid(); rng++) {
    // ideal mode just uses straigt line so we only
    // need 2 points per spline
    m_idealDACs[rng].resize(2);
    m_idealADCs[rng].resize(2);
    
    m_idealADCs[rng][0] = 0;
    m_idealADCs[rng][1] = maxADC;

    m_idealDACs[rng][0] = 0;
    m_idealDACs[rng][1] = 
      (unsigned int)(maxADC / owner->m_idealCalib.inlADCPerDAC[rng]);
  }

  // we don't have this info at this point
  // so why fake it?
  m_idealErr = 0;

  return StatusCode::SUCCESS;
}
