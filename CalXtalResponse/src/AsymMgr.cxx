// LOCAL INCLUDES
#include "AsymMgr.h"

// GLAST INCLUDES
#include "CalibData/Cal/Xpos.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// EXTLIB INCLUDES

// STD
#include <algorithm>
#include <typeinfo>

using namespace std;
using namespace CalDefs;
using namespace idents;

AsymMgr::AsymMgr(const IdealCalCalib &idealCalib) : 
  CalibItemMgr(CalibData::CAL_Asym, 
               idealCalib,
               N_SPLINE_TYPES) {

  // set size of spline lists (1 per xtal)
  for (unsigned i = 0; i < m_splineLists.size(); i++)
    m_splineLists[i].resize(XtalIdx::N_VALS);
};

bool AsymMgr::validateRangeBase(const CalXtalId&, CalibData::RangeBase *rngBase) {
  const vector<CalibData::ValSig> *asymLrg;
  const vector<CalibData::ValSig> *asymSm;
  const vector<CalibData::ValSig> *asymNSPB;
  const vector<CalibData::ValSig> *asymPSNB;

  CalibData::CalAsym *asym = (CalibData::CalAsym*)(rngBase);

  if (!(asymLrg = asym->getBig())) {
    // create MsgStream only when needed for performance
    MsgStream msglog(m_msgSvc, *m_logName); 
    msglog << MSG::ERROR << "Unable to retrieve calib data for " << m_calibPath << endreq;
    return false;
  }
  if (!(asymSm = asym->getSmall())) {
    // create MsgStream only when needed for performance
    MsgStream msglog(m_msgSvc, *m_logName); 
    msglog << MSG::ERROR << "Unable to retrieve calib data for " << m_calibPath << endreq;
    return false;
  }
  if (!(asymNSPB = asym->getNSmallPBig())) {
    // create MsgStream only when needed for performance
    MsgStream msglog(m_msgSvc, *m_logName); 
    msglog << MSG::ERROR << "Unable to retrieve calib data for " << m_calibPath << endreq;
    return false;
  }
  if (!(asymPSNB = asym->getPSmallNBig())) {
    // create MsgStream only when needed for performance
    MsgStream msglog(m_msgSvc, *m_logName); 
    msglog << MSG::ERROR << "Unable to retrieve calib data for " << m_calibPath << endreq;
    return false;
  }

  // get Xpos vals.
  unsigned XposSize= m_calibBase->getXpos()->getVals()->size();
  if (XposSize != asymLrg->size() ||
      XposSize != asymSm->size() ||
      XposSize != asymNSPB->size() ||
      XposSize != asymPSNB->size()) {
    // create MsgStream only when needed for performance
    MsgStream msglog(m_msgSvc, *m_logName); 
    msglog << MSG::ERROR << "Invalid # of vals for " << m_calibPath << endreq;
    return false;
  }
  return true;
}

/// retrieve Asymmetry calibration information for one xtal
StatusCode AsymMgr::getAsym(const CalXtalId &xtalId,
                            const vector<CalibData::ValSig> *&asymLrg,
                            const vector<CalibData::ValSig> *&asymSm,
                            const vector<CalibData::ValSig> *&asymNSPB,
                            const vector<CalibData::ValSig> *&asymPSNB,
                            const vector<float>  *&xVals) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  if (m_idealMode) {
    // return default vals if we're in ideal (fake) mode
    asymLrg  = &m_idealAsymLrg;
    asymSm   = &m_idealAsymSm;
    asymNSPB = &m_idealAsymNSPB;
    asymPSNB = &m_idealAsymPSNB;
    xVals    = &m_idealXVals;
    return StatusCode::SUCCESS;
  }
  
  CalibData::CalAsym *asym = getRangeBase(xtalId);
  
  // get main data arrays
  asymLrg = asym->getBig();
  asymSm = asym->getSmall();
  asymNSPB = asym->getNSmallPBig();
  asymPSNB = asym->getPSmallNBig();
  CalibData::Xpos *tmpXpos = m_calibBase->getXpos();
  xVals = tmpXpos->getVals();

  return StatusCode::SUCCESS;
}

StatusCode AsymMgr::genSplines() {
  StatusCode sc;

  // vector<double> arrays for input into genSpline
  vector<double> dblAsymLrg;
  vector<double> dblAsymSm;
  vector<double> dblAsymNSPB;
  vector<double> dblAsymPSNB;
  vector<double> dblXpos;

  for (XtalIdx xtalIdx; xtalIdx.isValid(); xtalIdx++) {
    // needed params for getAsym
    const vector<CalibData::ValSig> *asymLrg;
    const vector<CalibData::ValSig> *asymSm;
    const vector<CalibData::ValSig> *asymNSPB;
    const vector<CalibData::ValSig> *asymPSNB;
    const vector<float> *Xpos;

    sc = getAsym(xtalIdx.getCalXtalId(), 
                 asymLrg, asymSm, 
                 asymNSPB, asymPSNB, 
                 Xpos);
    if (sc.isFailure()) return sc;

    int n = Xpos->size();
    
    //////////////////////////////////////////
    //-- LINEARLY EXTRAPOLATED END-POINTS --//
    //////////////////////////////////////////

    // add one point to each end to ensure better
    // behaved spline functions
    
    dblAsymLrg.resize(n+2);
    dblAsymSm.resize(n+2);
    dblAsymPSNB.resize(n+2);
    dblAsymNSPB.resize(n+2);
    dblXpos.resize(n+2);

    for (int i = 0; i < n; i++) {
      dblAsymLrg[i+1]  = (*asymLrg)[i].getVal();
      dblAsymSm[i+1]   = (*asymSm)[i].getVal();
      dblAsymNSPB[i+1] = (*asymNSPB)[i].getVal();
      dblAsymPSNB[i+1] = (*asymPSNB)[i].getVal();

      dblXpos[i+1]     = (*Xpos)[i];
    }

    //-- LINEAR EXTRAPOLATION --//
    dblAsymLrg[0]  = 2*dblAsymLrg[1]  - dblAsymLrg[2];
    dblAsymSm[0]   = 2*dblAsymSm[1]   - dblAsymSm[2];
    dblAsymNSPB[0] = 2*dblAsymNSPB[1] - dblAsymNSPB[2];
    dblAsymPSNB[0] = 2*dblAsymPSNB[1] - dblAsymPSNB[2];

    dblXpos[0]     = 2*dblXpos[1] - dblXpos[2];

    dblAsymLrg[n+1]  = 2*dblAsymLrg[n]  - dblAsymLrg[n-1];
    dblAsymSm[n+1]   = 2*dblAsymSm[n]   - dblAsymSm[n-1];
    dblAsymNSPB[n+1] = 2*dblAsymNSPB[n] - dblAsymNSPB[n-1];
    dblAsymPSNB[n+1] = 2*dblAsymPSNB[n] - dblAsymPSNB[n-1];

    dblXpos[n+1]     = 2*dblXpos[n] - dblXpos[n-1];

    genSpline(ASYMLRG_SPLINE,   xtalIdx, "asymLrg",   dblXpos, dblAsymLrg);
    genSpline(ASYMSM_SPLINE,    xtalIdx, "asymSm",    dblXpos, dblAsymSm);
    genSpline(ASYMNSPB_SPLINE,  xtalIdx, "asymNSPB",  dblXpos, dblAsymNSPB);
    genSpline(ASYMPSNB_SPLINE,  xtalIdx, "asymPSNB",  dblXpos, dblAsymPSNB);

    genSpline(INV_ASYMLRG_SPLINE,   xtalIdx, "invLrg",   dblAsymLrg,   dblXpos);
    genSpline(INV_ASYMSM_SPLINE,    xtalIdx, "invSm",    dblAsymSm,    dblXpos);
    genSpline(INV_ASYMNSPB_SPLINE,  xtalIdx, "invNSPB",  dblAsymNSPB,  dblXpos);
    genSpline(INV_ASYMPSNB_SPLINE,  xtalIdx, "invPSNB",  dblAsymPSNB,  dblXpos);
  }  
  
  return StatusCode::SUCCESS;
}

StatusCode AsymMgr::fillRangeBases() {
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

StatusCode AsymMgr::loadIdealVals() {
  MsgStream msglog(m_msgSvc, *m_logName); 
  // linear 'fake' spline needs only 2 points
  m_idealAsymLrg.resize(2);
  m_idealAsymSm.resize(2);
  m_idealAsymNSPB.resize(2);
  m_idealAsymPSNB.resize(2);
  m_idealXVals.resize(2);

  m_idealAsymLrg[0].m_val = m_idealCalib.asymLrgNeg;
  m_idealAsymLrg[0].m_sig = m_idealCalib.asymLrgNeg * 
    m_idealCalib.asymSigPct;
  m_idealAsymLrg[1].m_val = m_idealCalib.asymLrgPos;
  m_idealAsymLrg[1].m_sig = m_idealCalib.asymLrgPos * 
    m_idealCalib.asymSigPct;

  m_idealAsymSm[0].m_val = m_idealCalib.asymSmNeg;
  m_idealAsymSm[0].m_sig = m_idealCalib.asymSmNeg * 
    m_idealCalib.asymSigPct;
  m_idealAsymSm[1].m_val = m_idealCalib.asymSmPos;
  m_idealAsymSm[1].m_sig = m_idealCalib.asymSmPos * 
    m_idealCalib.asymSigPct;

  m_idealAsymPSNB[0].m_val = m_idealCalib.asymPSNBNeg;
  m_idealAsymPSNB[0].m_sig = m_idealCalib.asymPSNBNeg * 
    m_idealCalib.asymSigPct;
  m_idealAsymPSNB[1].m_val = m_idealCalib.asymPSNBPos;
  m_idealAsymPSNB[1].m_sig = m_idealCalib.asymPSNBPos * 
    m_idealCalib.asymSigPct;

  m_idealAsymNSPB[0].m_val = m_idealCalib.asymNSPBNeg;
  m_idealAsymNSPB[0].m_sig = m_idealCalib.asymNSPBNeg * 
    m_idealCalib.asymSigPct;
  m_idealAsymNSPB[1].m_val = m_idealCalib.asymNSPBPos;
  m_idealAsymNSPB[1].m_sig = m_idealCalib.asymNSPBPos * 
    m_idealCalib.asymSigPct;

  float csiLength = 326.0;
  
  // 0 is at xtal center
  m_idealXVals[0] = -1*csiLength/2.0;
  m_idealXVals[1] = csiLength/2.0;
  
  return StatusCode::SUCCESS;
}
