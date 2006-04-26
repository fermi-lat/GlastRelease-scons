// $Header$
/** @file
    @author Zach Fewtrell
*/
// LOCAL INCLUDES
#include "AsymMgr.h"

// GLAST INCLUDES
#include "CalibData/Cal/Xpos.h"

// EXTLIB INCLUDES

// STD

using namespace std;
using namespace CalUtil;
using namespace idents;
using namespace CalibData;

/// used to represent invalid values in internal arrays.
const float BAD_FLOAT = -99999.9F;

AsymMgr::AsymMgr(CalCalibShared &ccsShared) : 
  CalibItemMgr(CAL_Asym,
               ccsShared,
               N_SPLINE_TYPES) 
{
  // set size of spline lists (1 per xtal)
  for (unsigned i = 0; i < m_splineLists.size(); i++) {
    m_splineLists[i].resize(XtalIdx::N_VALS, 0);
    m_splineXMin[i].resize(XtalIdx::N_VALS,  0);
    m_splineXMax[i].resize(XtalIdx::N_VALS,  0);
  }

  /// initialize asymCtr array
  for (AsymType asymType; asymType.isValid(); asymType++)
    m_asymCtr[asymType].fill(BAD_FLOAT);
}

/// get Asymmetry calibration information for one xtal
CalAsym *AsymMgr::getAsym(XtalIdx xtalIdx) {
  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return NULL;
  
  return (CalAsym*)m_rngBases[xtalIdx];
}

Xpos *AsymMgr::getXpos() {
  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return NULL;

  if (m_idealMode) return m_idealXpos.get();
  return m_calibBase->getXpos();
}

/** return p3 such that p3 - p2 = p2 - p1
 */
template <class Ty>
inline Ty extrap(Ty p1, Ty p2) {
  return 2*p2 - p1;
}

StatusCode AsymMgr::genLocalStore() {
  m_rngBases.resize(XtalIdx::N_VALS,0);

  StatusCode sc;

  // pointer to asym arrays in calib db
  CalArray<AsymType, const vector<ValSig> *> db_asym;
  // vector<float> arrays for input into each genSpline
  CalArray<AsymType, vector<float> > tmp_asym;
  vector<float> tmp_xvals;

  // used for description strings
  CalArray<AsymType, string> asymSuffix;
  asymSuffix[ASYM_LL] = "LL";
  asymSuffix[ASYM_LS] = "LS";
  asymSuffix[ASYM_SL] = "SL";
  asymSuffix[ASYM_SS] = "SS";

  for (XtalIdx xtalIdx; xtalIdx.isValid(); xtalIdx++) {
    if (!m_idealMode) {
      // RETRIEVE & VALIDATE RangeBase object from TDS
      CalAsym *asym = (CalAsym*)getRangeBase(xtalIdx.getCalXtalId());
      if (!validateRangeBase(asym)) continue;
      m_rngBases[xtalIdx] = asym;
    } else m_rngBases[xtalIdx] = m_idealAsym.get();

    // support missing towers & missing crystals
    // keep moving if we're missing a particular calibration
    const CalAsym *asym = getAsym(xtalIdx);
    if (!asym) continue; //support partial LATs

    db_asym[ASYM_LL] = asym->getBig();
    db_asym[ASYM_SS] = asym->getSmall();
    db_asym[ASYM_LS] = asym->getNSmallPBig();
    db_asym[ASYM_SL] = asym->getPSmallNBig();

    const Xpos *xpos = getXpos();
    if (!xpos) return StatusCode::FAILURE; // shouldn't happen
    const vector<float> *xvals = xpos->getVals();
    int n = xvals->size();
    
    //////////////////////////////////////////
    //-- LINEARLY EXTRAPOLATED END-POINTS --//
    //////////////////////////////////////////
    //-- LINEAR EXTRAPOLATION x points--//
    tmp_xvals.resize   (n + 4);
    for (int i = 0; i < n; i++) 
      tmp_xvals[i + 2] = (*xvals)[i];
      
    tmp_xvals[1]   = extrap(tmp_xvals[3],   tmp_xvals[2]);
    tmp_xvals[0]   = extrap(tmp_xvals[2],   tmp_xvals[1]);
    tmp_xvals[n+2] = extrap(tmp_xvals[n],   tmp_xvals[n+1]);
    tmp_xvals[n+3] = extrap(tmp_xvals[n+1], tmp_xvals[n+2]);

    for (AsymType asymType; asymType.isValid(); asymType++) {
      // add two points to each end to ensure better
      // good spline behavior past the edge of the
      // xtal.  if you think this isn't such a good
      // idea, then think about it for a while
      // or ask zach to explain.
      // 4 new points total
      
      tmp_asym[asymType].resize(n + 4);
      
      // load original points into middle of new, 
      // wider vector
      for (int i = 0; i < n; i++)
        tmp_asym[asymType][i+2] = (*db_asym[asymType])[i].getVal();

      // point 1
      tmp_asym[asymType][1] = extrap(tmp_asym[asymType][3], tmp_asym[asymType][2]);

      // point 0
      tmp_asym[asymType][0] = extrap(tmp_asym[asymType][2], tmp_asym[asymType][1]);

      // 2nd last point
      tmp_asym[asymType][n+2] = extrap(tmp_asym[asymType][n], tmp_asym[asymType][n+1]);
    
      // last point
      tmp_asym[asymType][n+3] = extrap(tmp_asym[asymType][n+1], tmp_asym[asymType][n+2]);

      // put xtal id string into spline name
      ostringstream xtalStr;
      xtalStr << '[' << xtalIdx.getCalXtalId() << ']';
      
      string splineName = "asym" + asymSuffix[asymType];
      genSpline(splineIdx(asymType, POS2ASYM),  xtalIdx, splineName  + xtalStr.str(),   
                tmp_xvals, tmp_asym[asymType]);

      splineName = "invAsym" + asymSuffix[asymType];
      genSpline(splineIdx(asymType, ASYM2POS),  xtalIdx, splineName + xtalStr.str(),   
                tmp_asym[asymType],   tmp_xvals);

      //-- ASYM CTR --//
      float asymCtr;
      sc = evalAsym(xtalIdx, asymType, 0, asymCtr);
      if (sc.isFailure()) return sc;
      m_asymCtr[asymType][xtalIdx] = asymCtr;
    }  
  }  

  return StatusCode::SUCCESS;
}

StatusCode AsymMgr::loadIdealVals() {
  // linear 'fake' spline needs only 2 points
  vector<ValSig> big(2);
  vector<ValSig> small(2);
  vector<ValSig> nspb(2);
  vector<ValSig> psnb(2);

  big[0].m_val = m_ccsShared.m_idealCalib.asymLrgNeg;
  big[0].m_sig = m_ccsShared.m_idealCalib.asymLrgNeg * 
    m_ccsShared.m_idealCalib.asymSigPct;
  big[1].m_val = m_ccsShared.m_idealCalib.asymLrgPos;
  big[1].m_sig = m_ccsShared.m_idealCalib.asymLrgPos * 
    m_ccsShared.m_idealCalib.asymSigPct;

  small[0].m_val = m_ccsShared.m_idealCalib.asymSmNeg;
  small[0].m_sig = m_ccsShared.m_idealCalib.asymSmNeg * 
    m_ccsShared.m_idealCalib.asymSigPct;
  small[1].m_val = m_ccsShared.m_idealCalib.asymSmPos;
  small[1].m_sig = m_ccsShared.m_idealCalib.asymSmPos * 
    m_ccsShared.m_idealCalib.asymSigPct;

  psnb[0].m_val = m_ccsShared.m_idealCalib.asymPSNBNeg;
  psnb[0].m_sig = m_ccsShared.m_idealCalib.asymPSNBNeg * 
    m_ccsShared.m_idealCalib.asymSigPct;
  psnb[1].m_val = m_ccsShared.m_idealCalib.asymPSNBPos;
  psnb[1].m_sig = m_ccsShared.m_idealCalib.asymPSNBPos * 
    m_ccsShared.m_idealCalib.asymSigPct;

  nspb[0].m_val = m_ccsShared.m_idealCalib.asymNSPBNeg;
  nspb[0].m_sig = m_ccsShared.m_idealCalib.asymNSPBNeg * 
    m_ccsShared.m_idealCalib.asymSigPct;
  nspb[1].m_val = m_ccsShared.m_idealCalib.asymNSPBPos;
  nspb[1].m_sig = m_ccsShared.m_idealCalib.asymNSPBPos * 
    m_ccsShared.m_idealCalib.asymSigPct;

  m_idealAsym.reset(new CalAsym(&big, &small, &nspb, &psnb));

  
  // 0 is at xtal center
  vector<float> xvals(2);
  
  xvals[0] = -1*m_csiLength/2.0;
  xvals[1] = m_csiLength/2.0;

  m_idealXpos.reset(new Xpos(&xvals));
  
  return StatusCode::SUCCESS;
}

bool AsymMgr::validateRangeBase(CalAsym *asym) {
  if (!asym) return false;

  const vector<ValSig> *asym_ll;
  const vector<ValSig> *asym_ss;
  const vector<ValSig> *asym_ls;
  const vector<ValSig> *asym_sl;

  if (!(asym_ll = asym->getBig())) {
    // no error print out req'd b/c we're supporting LAT configs w/ empty bays
    // however, if asym->getBig() is successful & following checks fail
    // then we have a problem b/c we have calib data which is only good for
    // partial xtal.
    return false;
  }
  if (!(asym_ss = asym->getSmall())        ||
      !(asym_ls = asym->getNSmallPBig()) ||
      !(asym_sl = asym->getPSmallNBig())) {
    // create MsgStream only when needed for performance
    MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name()); 
    msglog << MSG::ERROR << "can't get calib data for " 
           << m_calibPath;
    msglog << endreq;
    return false;
  }

  // get Xpos vals.
  unsigned XposSize= getXpos()->getVals()->size();
  if (XposSize != asym_ll->size() ||
      XposSize != asym_ss->size() ||
      XposSize != asym_ls->size() ||
      XposSize != asym_sl->size()) {
    // create MsgStream only when needed for performance
    MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name()); 
    msglog << MSG::ERROR << "Invalid # of vals for " << m_calibPath << endreq;
    return false;
  }
  return true;
}

StatusCode AsymMgr::initialize(const string &flavor) {
  StatusCode sc;

  // try to find the GlastDevSvc service
  IGlastDetSvc* detSvc;
  sc = m_ccsShared.m_service->service("GlastDetSvc", detSvc);
  if (sc.isFailure() ) {
    MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name()); 
    msglog << MSG::ERROR << "  can't get GlastDetSvc " << endreq;
    return sc;
  }

  double tmp;
  sc = detSvc->getNumericConstByName("CsILength", &tmp);
  m_csiLength = tmp;
  if (sc.isFailure()) {
    MsgStream msglog(m_ccsShared.m_service->msgSvc(), m_ccsShared.m_service->name()); 
    msglog << MSG::ERROR << " constant CsILength not defined" << endreq;
    return StatusCode::FAILURE;
  }

  /// call chain back to parent class routine as well
  return CalibItemMgr::initialize(flavor);
}


StatusCode AsymMgr::getAsymCtr(XtalIdx xtalIdx, AsymType asymType, float &asymCtr) {
  // make sure we have valid calib data for this event.
  StatusCode sc;
  sc = updateCalib();
  if (sc.isFailure()) return StatusCode::FAILURE;

  asymCtr = m_asymCtr[asymType][xtalIdx];

  return StatusCode::SUCCESS;
}
