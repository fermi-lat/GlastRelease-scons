// LOCAL
#include "IntNonlinMgr.h"

// GLAST
#include "CalibData/DacCol.h"

// EXTLIB

// STD
#include <algorithm>

using namespace CalDefs;
using namespace idents;

IntNonlinMgr::IntNonlinMgr(const IdealCalCalib &idealCalib) : 
  CalibItemMgr(CalibData::CAL_IntNonlin, 
               idealCalib, 
               N_SPLINE_TYPES),
  m_idealADCs(RngNum::N_VALS), // one spline per range
  m_idealDACs(RngNum::N_VALS)  // one spline per range
{

  // set size of spline lists (1 per range)
  for (unsigned i = 0; i < m_splineLists.size(); i++)
    m_splineLists[i].resize(RngIdx::N_VALS);
};

bool IntNonlinMgr::validateRangeBase(const CalXtalId &xtalId, CalibData::RangeBase *rngBase) {
  // recast for specific calibration type
  CalibData::IntNonlin* intNonlin = (CalibData::IntNonlin*)(rngBase);

  //get vector of vals
  const vector<float> *intNonlinVec = intNonlin->getValues();
  if (!intNonlinVec) {
    // create MsgStream only when needed for performance
    MsgStream msglog(m_msgSvc, *m_logName); 
    msglog << MSG::ERROR << "Unable to get vector for IntNonlin vals" << endreq;
    return false;
  }

  //get collection of associated DAC vals
  CalibData::DacCol *intNonlinDacCol = m_calibBase->getDacCol(xtalId.getRange());
  const vector<unsigned> *intNonlinDacVec;
  if (intNonlinDacCol) {
    intNonlinDacVec = intNonlinDacCol->getDacs();
  } else {
    // create MsgStream only when needed for performance
    MsgStream msglog(m_msgSvc, *m_logName); 
    msglog << MSG::ERROR << "No intNonlinDacCol found.\n" << endreq;
    return false;
  }

  // check that there are enough DAC vals for ADC vals in this rng
  if (intNonlinVec->size() > intNonlinDacVec->size()) {
    // create MsgStream only when needed for performance
    MsgStream msglog(m_msgSvc, *m_logName); 
    msglog << MSG::ERROR << "intNonlin nADC > nDAC" 
           << " CalXtalId=" << xtalId
           << " nADC="   << intNonlinVec->size()   << " " 
           << " nDAC="   << intNonlinDacVec->size() << "CalXtalId=" << xtalId
           << endreq;
    return false;
  }

  return true;
}

/// retrieve integral non-linearity vals for given xtal/face/rng
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

  CalibData::IntNonlin *intNonlin = getRangeBase(xtalId);
  if (!intNonlin) return StatusCode::FAILURE;

  //get vector of vals
  const vector<float> *intNonlinVec = intNonlin->getValues();

  //get collection of associated DAC vals
  CalibData::DacCol *intNonlinDacCol = m_calibBase->getDacCol(xtalId.getRange());
  const vector<unsigned> *intNonlinDacVec;
  intNonlinDacVec = intNonlinDacCol->getDacs();

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

    sc = getIntNonlin(rngIdx.getCalXtalId(), adc, dac, error);
    if (sc.isFailure()) return sc;

    int n = min(adc->size(),dac->size());

    dblDac.resize(n);
    dblAdc.resize(n);

    // create double vector for input to genSpline()
    copy(adc->begin(), adc->begin() + n, dblAdc.begin());
    copy(dac->begin(), dac->begin() + n, dblDac.begin());

    genSpline(INL_SPLINE,     rngIdx, "intNonlin",    dblAdc, dblDac);
    genSpline(INV_INL_SPLINE, rngIdx, "invIntNonlin", dblDac, dblAdc);
  }

  return StatusCode::SUCCESS;
}

StatusCode IntNonlinMgr::fillRangeBases() {
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

StatusCode IntNonlinMgr::loadIdealVals() {
  MsgStream msglog(m_msgSvc, *m_logName); 

  //-- SANITY CHECKS --//
  if (m_idealCalib.ciULD.size() != (unsigned)RngNum::N_VALS) {
    msglog << MSG::ERROR << "Bad # of ULD vals in ideal CalCalib xml file" 
           << endreq;
    return StatusCode::FAILURE;
  }
  if (m_idealCalib.inlADCPerDAC.size() != (unsigned)RngNum::N_VALS) {
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
    m_idealADCs[rng][1] = m_idealCalib.ciULD[rng];

    m_idealDACs[rng][0] = 0;
    m_idealDACs[rng][1] = 
      (unsigned int)(m_idealCalib.ciULD[rng] /
                     m_idealCalib.inlADCPerDAC[rng]);
  }

  // we don't have this info at this point
  // so why fake it?
  m_idealErr = 0;

  return StatusCode::SUCCESS;
}
