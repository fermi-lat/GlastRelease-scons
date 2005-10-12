// LOCAL INCLUDES
#include "AsymMgr.h"
#include "CalCalibSvc.h"

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

AsymMgr::AsymMgr() : 
  CalibItemMgr(CalibData::CAL_Asym, 
               N_SPLINE_TYPES) {

  // set size of spline lists (1 per xtal)
  for (unsigned i = 0; i < m_splineLists.size(); i++) {
    m_splineLists[i].resize(XtalIdx::N_VALS, 0);
    m_splineXMin[i].resize(XtalIdx::N_VALS,  0);
    m_splineXMax[i].resize(XtalIdx::N_VALS,  0);
  }
}

/// get Asymmetry calibration information for one xtal
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

  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return sc;

  
  CalibData::CalAsym *asym = 
	  (CalibData::CalAsym *)getRangeBase(xtalId);
  if (!asym) return StatusCode::FAILURE;
  
  // get main data arrays
  asymLrg = asym->getBig();
  asymSm = asym->getSmall();
  asymNSPB = asym->getNSmallPBig();
  asymPSNB = asym->getPSmallNBig();
  CalibData::Xpos *tmpXpos = m_calibBase->getXpos();
  xVals = tmpXpos->getVals();

  return StatusCode::SUCCESS;
}

/** return p3 such that p3 - p2 = p2 - p1
 */
inline double extrap(double p1, double p2) {
  return 2*p2 - p1;
}

StatusCode AsymMgr::genLocalStore() {
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

    // support missing towers & missing crystals
    // keep moving if we're missing a particular calibration
    sc = getAsym(xtalIdx.getCalXtalId(), 
                 asymLrg, asymSm, 
                 asymNSPB, asymPSNB, 
                 Xpos);
    if (sc.isFailure()) continue; //support partial LATs

    int n = Xpos->size();
    
    //////////////////////////////////////////
    //-- LINEARLY EXTRAPOLATED END-POINTS --//
    //////////////////////////////////////////

    // add two points to each end to ensure better
    // good spline behavior past the edge of the
    // xtal.  if you think this isn't such a good
    // idea, then think about it for a while
    // or ask zach to explain.
    // 4 new points total
    
    dblAsymLrg.resize (n + 4);
    dblAsymSm.resize  (n + 4);
    dblAsymPSNB.resize(n + 4);
    dblAsymNSPB.resize(n + 4);
    dblXpos.resize    (n + 4);

    // load original points into middle of new, 
    // wider vector
    for (int i = 0; i < n; i++) {
      dblAsymLrg [i + 2] = (*asymLrg) [i].getVal();
      dblAsymSm  [i + 2] = (*asymSm)  [i].getVal();
      dblAsymNSPB[i + 2] = (*asymNSPB)[i].getVal();
      dblAsymPSNB[i + 2] = (*asymPSNB)[i].getVal();

      dblXpos    [i + 2] = (*Xpos)[i];
    }

    //-- LINEAR EXTRAPOLATION --//
    // point 1
    dblAsymLrg [1] = extrap(dblAsymLrg [3], dblAsymLrg [2]);
    dblAsymSm  [1] = extrap(dblAsymSm  [3], dblAsymSm  [2]);
    dblAsymNSPB[1] = extrap(dblAsymNSPB[3], dblAsymNSPB[2]);
    dblAsymPSNB[1] = extrap(dblAsymPSNB[3], dblAsymPSNB[2]);
    dblXpos    [1] = extrap(dblXpos    [3], dblXpos    [2]);

    // point 0
    dblAsymLrg [0] = extrap(dblAsymLrg [2], dblAsymLrg [1]);
    dblAsymSm  [0] = extrap(dblAsymSm  [2], dblAsymSm  [1]);
    dblAsymNSPB[0] = extrap(dblAsymNSPB[2], dblAsymNSPB[1]);
    dblAsymPSNB[0] = extrap(dblAsymPSNB[2], dblAsymPSNB[1]);
    dblXpos    [0] = extrap(dblXpos    [2], dblXpos    [1]);

    // 2nd last point
    dblAsymLrg [n+2] = extrap(dblAsymLrg [n], dblAsymLrg [n+1]);
    dblAsymSm  [n+2] = extrap(dblAsymSm  [n], dblAsymSm  [n+1]);
    dblAsymNSPB[n+2] = extrap(dblAsymNSPB[n], dblAsymNSPB[n+1]);
    dblAsymPSNB[n+2] = extrap(dblAsymPSNB[n], dblAsymPSNB[n+1]);
    dblXpos    [n+2] = extrap(dblXpos    [n], dblXpos    [n+1]);
    
    // last point
    dblAsymLrg [n+3] = extrap(dblAsymLrg [n+1], dblAsymLrg [n+2]);
    dblAsymSm  [n+3] = extrap(dblAsymSm  [n+1], dblAsymSm  [n+2]);
    dblAsymNSPB[n+3] = extrap(dblAsymNSPB[n+1], dblAsymNSPB[n+2]);
    dblAsymPSNB[n+3] = extrap(dblAsymPSNB[n+1], dblAsymPSNB[n+2]);
    dblXpos    [n+3] = extrap(dblXpos    [n+1], dblXpos    [n+2]);

    if (owner->m_superVerbose) {
      // create MsgStream only when needed (for performance)
      MsgStream msglog(owner->msgSvc(), owner->name()); 
      msglog << MSG::VERBOSE << "xpos ";
      for (unsigned i = 0; i < dblXpos.size(); i++) 
        msglog << dblXpos[i] << ' ';
      msglog << endreq;

      msglog << MSG::VERBOSE << "asymLL ";
      for (unsigned i = 0; i < dblAsymLrg.size(); i++) 
        msglog << dblAsymLrg[i] << ' ';
      msglog << endreq;
    }

	CalXtalId xtalId = xtalIdx.getCalXtalId();
	
	// put xtal id string into spline name
    ostringstream xtalStr;
    xtalStr << '[' << xtalId << ']';

    genSpline(ASYMLRG_SPLINE,   xtalId, "asymLrg"     + xtalStr.str(),   
              dblXpos, dblAsymLrg);
    genSpline(ASYMSM_SPLINE,    xtalId, "asymSm"      + xtalStr.str(),    
              dblXpos, dblAsymSm);
    genSpline(ASYMNSPB_SPLINE,  xtalId, "asymNSPB"    + xtalStr.str(),  
              dblXpos, dblAsymNSPB);
    genSpline(ASYMPSNB_SPLINE,  xtalId, "asymPSNB"    + xtalStr.str(),  
              dblXpos, dblAsymPSNB);

    genSpline(INV_ASYMLRG_SPLINE,   xtalId, "invAsymLrg"  + xtalStr.str(),   
              dblAsymLrg,   dblXpos);
    genSpline(INV_ASYMSM_SPLINE,    xtalId, "invAsymSm"   + xtalStr.str(),    
              dblAsymSm,    dblXpos);
    genSpline(INV_ASYMNSPB_SPLINE,  xtalId, "invAsymNSPB" + xtalStr.str(),  
              dblAsymNSPB,  dblXpos);
    genSpline(INV_ASYMPSNB_SPLINE,  xtalId, "invAsymPSNB" + xtalStr.str(),  
              dblAsymPSNB,  dblXpos);
  }  
  
  return StatusCode::SUCCESS;
}

StatusCode AsymMgr::loadIdealVals() {
  // linear 'fake' spline needs only 2 points
  m_idealAsymLrg.resize (2);
  m_idealAsymSm.resize  (2);
  m_idealAsymNSPB.resize(2);
  m_idealAsymPSNB.resize(2);
  m_idealXVals.resize   (2);

  m_idealAsymLrg[0].m_val = owner->m_idealCalib.asymLrgNeg;
  m_idealAsymLrg[0].m_sig = owner->m_idealCalib.asymLrgNeg * 
    owner->m_idealCalib.asymSigPct;
  m_idealAsymLrg[1].m_val = owner->m_idealCalib.asymLrgPos;
  m_idealAsymLrg[1].m_sig = owner->m_idealCalib.asymLrgPos * 
    owner->m_idealCalib.asymSigPct;

  m_idealAsymSm[0].m_val = owner->m_idealCalib.asymSmNeg;
  m_idealAsymSm[0].m_sig = owner->m_idealCalib.asymSmNeg * 
    owner->m_idealCalib.asymSigPct;
  m_idealAsymSm[1].m_val = owner->m_idealCalib.asymSmPos;
  m_idealAsymSm[1].m_sig = owner->m_idealCalib.asymSmPos * 
    owner->m_idealCalib.asymSigPct;

  m_idealAsymPSNB[0].m_val = owner->m_idealCalib.asymPSNBNeg;
  m_idealAsymPSNB[0].m_sig = owner->m_idealCalib.asymPSNBNeg * 
    owner->m_idealCalib.asymSigPct;
  m_idealAsymPSNB[1].m_val = owner->m_idealCalib.asymPSNBPos;
  m_idealAsymPSNB[1].m_sig = owner->m_idealCalib.asymPSNBPos * 
    owner->m_idealCalib.asymSigPct;

  m_idealAsymNSPB[0].m_val = owner->m_idealCalib.asymNSPBNeg;
  m_idealAsymNSPB[0].m_sig = owner->m_idealCalib.asymNSPBNeg * 
    owner->m_idealCalib.asymSigPct;
  m_idealAsymNSPB[1].m_val = owner->m_idealCalib.asymNSPBPos;
  m_idealAsymNSPB[1].m_sig = owner->m_idealCalib.asymNSPBPos * 
    owner->m_idealCalib.asymSigPct;

  float csiLength = 326.0;
  
  // 0 is at xtal center
  m_idealXVals[0] = -1*csiLength/2.0;
  m_idealXVals[1] = csiLength/2.0;
  
  return StatusCode::SUCCESS;
}
