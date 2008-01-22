// $Header$
/** @file
    @author Z.Fewtrell
*/
// LOCAL
#include "IntNonlinMgr.h"

// GLAST
#include "CalibData/DacCol.h"

// EXTLIB

// STD
#include <algorithm>

using namespace CalUtil;
using namespace idents;

IntNonlinMgr::IntNonlinMgr(CalCalibShared &ccsShared) : 
  CalibItemMgr(ICalibPathSvc::Calib_CAL_IntNonlin, 
               ccsShared,
               CalUtil::RngIdx::N_VALS,
               N_SPLINE_TYPES)
{
}

const vector<float> *IntNonlinMgr::getInlAdc(const CalUtil::RngIdx rngIdx) {
  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return NULL;

  CalibData::IntNonlin const*const inl = (CalibData::IntNonlin*)m_rngBases[rngIdx];

  if (!inl) return NULL;

  return inl->getValues();

}

const vector<float> *IntNonlinMgr::getInlCIDAC(const CalUtil::RngIdx rngIdx) {
  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return NULL;

  // return 0 if array is empty
  if (m_CIDACs[rngIdx].size() == 0)
	  return NULL;

  return &(m_CIDACs[rngIdx]);
}


/** return y3 such that (y2 - y1)/(x2 - x1) = (y3 - y2)/(x3 - x2)
 */
template <class Ty>
inline Ty extrap(Ty x1, Ty x2, Ty x3, Ty y1, Ty y2) {
  return (x3-x2)*(y2-y1)/(x2-x1) + y2;
}


StatusCode IntNonlinMgr::genLocalStore() {
  for (RngIdx rngIdx; rngIdx.isValid(); rngIdx++) {
    const vector<float> *adc;
    const RngNum rng = rngIdx.getRng();

    //-- IDEAL MODE --//
    if (m_idealMode) {
      adc = m_idealINL[rng].get()->getValues();
      const vector<float> *cidac = m_idealINL[rng].get()->getSdacs();

      m_CIDACs[rngIdx].resize(cidac->size());
      copy(cidac->begin(),
           cidac->end(),
           m_CIDACs[rngIdx].begin());

	  m_rngBases[rngIdx] = m_idealINL[rng].get();
    }
    
    //-- NORMAL (NON-IDEAL) MODE -//
    else {
      CalibData::IntNonlin const*const intNonlin 
        = (CalibData::IntNonlin *)getRangeBase(rngIdx.getCalXtalId());
      // support partial LAT
      if (!intNonlin) continue;
      if (!validateRangeBase(intNonlin)) continue;
      m_rngBases[rngIdx] = intNonlin;

      // quick check that ADC data is present 
      adc = intNonlin->getValues();
      // support partial LAT
      if (!adc) continue;
    
      //-- PHASE 1: populate CIDAC values (needed for spline generation )
      //-- 1st choice, use per-channel 'sdac' info if present
      const vector<float> *cidac = intNonlin->getSdacs();
      if (cidac) {
        m_CIDACs[rngIdx].resize(cidac->size());
        copy(cidac->begin(),
             cidac->end(),
             m_CIDACs[rngIdx].begin());
      }

      //-- 2nd choise, fall back to global 'DacCol' info
      else {
        //get collection of associated DAC vals
        CalibData::DacCol const*const intNonlinDacCol = 
          m_calibBase->getDacCol((CalXtalId::AdcRange)rng);
        
        const vector<unsigned> *globalCIDACs;
        globalCIDACs = intNonlinDacCol->getDacs();

        // if we've gotten this far, then we need 
        // the data to be present. can't have ADC
        // values for a channel w/ no matchin DAC
        // info
        if (!globalCIDACs) {
          // create MsgStream only when needed (for performance)
          MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name()); 
          msglog << MSG::ERROR << "ADC data w/ no matching CIDAC data (either per-channel or global) ADC channel: "
                 << rngIdx.val() 
                 << endreq;
          return StatusCode::FAILURE;
        }
        
        m_CIDACs[rngIdx].resize(globalCIDACs->size());
        copy(globalCIDACs->begin(),
             globalCIDACs->end(),
             m_CIDACs[rngIdx].begin());
      }
    }

    //-- check that we have enough points
    const vector<float> &cidac = m_CIDACs[rngIdx];
    if (cidac.size() < adc->size()) {
      MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name()); 
      msglog << MSG::ERROR << "Not enough CIDAC points for ADC channel: " << rngIdx.val()
             << endreq;
      return StatusCode::FAILURE;
    }

    //-- PHASE 2: generate splines
    vector<float> sp_cidac;
    vector<float> sp_adc;
    
    int n = min(adc->size(),cidac.size());

    // we need extra point for low end extrapolation
    sp_cidac.resize(n+1);
    sp_adc.resize(n+1);

    // create float vector for input to genSpline()
    // leave point at beginning of vector for low end extrapolation
    copy(adc->begin(),  adc->begin() + n,  sp_adc.begin()+1);
    copy(cidac.begin(), cidac.begin() + n, sp_cidac.begin()+1);

    // EXTRAPOLATE
    sp_cidac[0] = -200;
    sp_adc[0] = extrap(sp_cidac[2], sp_cidac[1], sp_cidac[0],
                       sp_adc[2], sp_adc[1]);

    // put rng id string into spline name
    ostringstream rngStr;
    rngStr << '[' << rngIdx.getCalXtalId()
           << ']';

    genSpline(INL_SPLINE, rngIdx, "INL" + rngStr.str(),    
              sp_adc, sp_cidac);
    genSpline(INV_INL_SPLINE, rngIdx, "invINL" + rngStr.str(), 
              sp_cidac, sp_adc);
  }

  return StatusCode::SUCCESS;
}

StatusCode IntNonlinMgr::loadIdealVals() {
  const unsigned short maxADC = 4095;

  //-- SANITY CHECKS --//
  if (m_ccsShared.m_idealCalib.inlADCPerCIDAC.size() != RngNum::N_VALS) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name()); 
    msglog << MSG::ERROR << "Bad # of ADCPerCIDAC vals in ideal CalCalib xml file" 
           << endreq;
    return StatusCode::FAILURE;
  }

  // ideal mode just uses straigt line so we only
  // need 2 points per spline
  vector<float> idealADCs(2);
  vector<float> idealCIDACs(2);
  for (RngNum rng; rng.isValid(); rng++) {
    // inl spline is pedestal subtracted.
    const float pedestal = m_ccsShared.m_idealCalib.pedVals[rng.val()];
    float maxADCPed = maxADC - pedestal;

    // optionally clip MAX ADC in HEX1 range to ULD value supplied by tholdci
    if (rng == HEX1)
      maxADCPed = min<float>(maxADCPed, m_ccsShared.m_idealCalib.ciULD[HEX1.val()]);

    idealADCs[0] = 0;
    idealADCs[1] = maxADCPed;

    idealCIDACs[0] = 0;
    idealCIDACs[1] = 
      maxADCPed / m_ccsShared.m_idealCalib.inlADCPerCIDAC[rng.val()];

    m_idealINL[rng].reset(new CalibData::IntNonlin(&idealADCs, 0, &idealCIDACs));
  }

  return StatusCode::SUCCESS;
}


bool IntNonlinMgr::validateRangeBase(CalibData::IntNonlin const*const intNonlin) {
  if (!intNonlin) return false;

  //get vector of vals
  const vector<float> *intNonlinVec = intNonlin->getValues();
  if (!intNonlinVec)    
    return false;

  return true;
}


StatusCode IntNonlinMgr::evalCIDAC(const CalUtil::RngIdx rngIdx, const float adc, float &cidac) {
    if (evalSpline(INL_SPLINE, rngIdx, adc, cidac).isFailure())
        return StatusCode::FAILURE;

    // ceiling check
    cidac = min(m_splineYMax[INL_SPLINE][rngIdx],cidac);

    return StatusCode::SUCCESS;
  }
  
StatusCode IntNonlinMgr::evalADC(const CalUtil::RngIdx rngIdx, const float cidac, float &adc) {
    if (evalSpline(INV_INL_SPLINE, rngIdx, cidac, adc).isFailure())
        return StatusCode::FAILURE;

    // ceiling check
    adc = min(m_splineYMax[INV_INL_SPLINE][rngIdx],adc);

    return StatusCode::SUCCESS;
}
