// $Header$
/** @file
    @author Zach Fewtrell
 */
// LOCAL
#include "IntNonlinMgr.h"
#include "CalCalibSvc.h"

// GLAST
#include "CalibData/DacCol.h"

// EXTLIB

// STD
#include <algorithm>

using namespace CalUtil;
using namespace idents;

IntNonlinMgr::IntNonlinMgr() : 
  CalibItemMgr(CalibData::CAL_IntNonlin, 
               N_SPLINE_TYPES),
  m_idealErr(0)
{
  // set size of spline lists (1 per range)
  for (unsigned i = 0; i < m_splineLists.size(); i++) {
    m_splineLists[i].resize(RngIdx::N_VALS, 0);
    m_splineXMin[i].resize (RngIdx::N_VALS, 0);
    m_splineXMax[i].resize (RngIdx::N_VALS, 0);
  }
}

/// get integral non-linearity vals for given xtal/face/rng
StatusCode IntNonlinMgr::getIntNonlin(RngIdx rngIdx,
                                      const vector< float > *&adcs,
                                      const vector< float > *&dacs,
                                      float &error) {
  if (m_idealMode) {
    RngNum rng(rngIdx.getRng());
    adcs = &m_idealADCs[rng];
    dacs = &m_idealDACs[rng];
    error = m_idealErr;
    return StatusCode::SUCCESS;
  }

  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return sc;

  
  CalibData::IntNonlin *intNonlin 
	  = (CalibData::IntNonlin*)m_rngBases[rngIdx];
  if (!intNonlin) return StatusCode::FAILURE;

  //-- retrieve ADC vals
  adcs = intNonlin->getValues();
  
  //-- retrieve DAC vals
  dacs = &(m_DACs[rngIdx]);

  // check array lens
  if (adcs->size() > dacs->size())
    return StatusCode::FAILURE;

  //whew, we made it this far... let's assign our outputs & leave!
  error = intNonlin->getError();

  return StatusCode::SUCCESS;
}

StatusCode IntNonlinMgr::genLocalStore() {
  const vector<float> *adc;
  vector<float> *dac;

  m_rngBases.resize(RngIdx::N_VALS, 0);

  for (RngIdx rngIdx; rngIdx.isValid(); rngIdx++) {
    RngNum rng = rngIdx.getRng();

    //-- IDEAL MODE --//
    if (m_idealMode) {
      adc = &m_idealADCs[rng];
      dac = &m_idealDACs[rng];
    }
    
    //-- NORMAL (NON-IDEAL) MODE -//
    else {
      dac = &m_DACs[rngIdx];
      
      CalibData::IntNonlin *intNonlin 
        = (CalibData::IntNonlin *)getRangeBase(rngIdx.getCalXtalId());
      // support partial LAT
      if (!intNonlin) continue;
      if (!validateRangeBase(intNonlin)) continue;
      m_rngBases[rngIdx] = intNonlin;

      // quick check that ADC data is present 
      adc = intNonlin->getValues();
      // support partial LAT
      if (!adc) continue;
    
      //-- PHASE 1: populate DAC values (needed for spline generation )
      //-- 1st choice, use per-channel 'sdac' info if present
      const vector<float> *sdacs = intNonlin->getSdacs();
      if (sdacs) {
        dac->resize(sdacs->size());
        copy(sdacs->begin(),
             sdacs->end(),
             dac->begin());
      }

      //-- 2nd choise, fall back to global 'DacCol' info
      else {
        //get collection of associated DAC vals
        CalibData::DacCol *intNonlinDacCol = 
          m_calibBase->getDacCol((CalXtalId::AdcRange)rng);
        
        const vector<unsigned> *globalDACs;
        globalDACs = intNonlinDacCol->getDacs();

        // if we've gotten this far, then we need 
        // the data to be present. can't have ADC
        // values for a channel w/ no matchin DAC
        // info
        if (!globalDACs) {
          // create MsgStream only when needed (for performance)
          MsgStream msglog(owner->msgSvc(), owner->name()); 
          msglog << MSG::ERROR << "ADC data w/ no matching DAC data (either per-channel or global"
                 << endreq;
          return StatusCode::FAILURE;
        }
        
        dac->resize(globalDACs->size());
        copy(globalDACs->begin(),
             globalDACs->end(),
             dac->begin());
      }
    }

    //-- PHASE 2: generate splines
    vector<float> dblDac;
    vector<float> dblAdc;
    
    int n = min(adc->size(),dac->size());

    dblDac.resize(n);
    dblAdc.resize(n);

    // create float vector for input to genSpline()
    copy(adc->begin(), adc->begin() + n, dblAdc.begin());
    copy(dac->begin(), dac->begin() + n, dblDac.begin());

    // put rng id string into spline name
    ostringstream rngStr;
    rngStr << '[' << rngIdx.getCalXtalId()
           << ']';

    genSpline(INL_SPLINE, rngIdx, "INL"    + rngStr.str(),    
              dblAdc, dblDac);
    genSpline(INV_INL_SPLINE, rngIdx, "invINL" + rngStr.str(), 
              dblDac, dblAdc);
  }

  return StatusCode::SUCCESS;
}

StatusCode IntNonlinMgr::loadIdealVals() {
  const int maxADC = 4095;

  //-- SANITY CHECKS --//
  if (owner->m_idealCalib.ciULD.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::ERROR << "Bad # of ULD vals in ideal CalCalib xml file" 
           << endreq;
    return StatusCode::FAILURE;
  }
  if (owner->m_idealCalib.inlADCPerDAC.size() != (unsigned)RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
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
      (unsigned int)(maxADC / owner->m_idealCalib.inlADCPerDAC[rng.getInt()]);
  }

  // we don't have this info at this point
  // so why fake it?
  m_idealErr = 0;

  return StatusCode::SUCCESS;
}


bool IntNonlinMgr::validateRangeBase(CalibData::IntNonlin *intNonlin) {
  if (!intNonlin) return false;

  //get vector of vals
  const vector<float> *intNonlinVec = intNonlin->getValues();
  if (!intNonlinVec)    
    return false;

  return true;
}
