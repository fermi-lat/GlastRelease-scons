// LOCAL INCLUDES
#include "AsymMgr.h"

// GLAST INCLUDES
#include "CalibData/Cal/Xpos.h"

// EXTLIB INCLUDES

// STD
#include <algorithm>
#include <typeinfo>

using namespace std;
using namespace CalDefs;
using namespace idents;

AsymMgr::AsymMgr() : 
  CalibItemMgr(CalibData::CAL_Asym, N_SPLINE_TYPES) {

  // set size of spline lists (1 per xtal)
  for (unsigned i = 0; i < m_splineLists.size(); i++)
    m_splineLists[i].resize(XtalIdx::N_VALS);
};

bool AsymMgr::validateRangeBase(const CalXtalId &xtalId, CalibData::RangeBase *rngBase) {
  MsgStream msglog(m_msgSvc, *m_logName);
  const vector<CalibData::ValSig> *lrg;
  const vector<CalibData::ValSig> *sm;
  const vector<CalibData::ValSig> *nspb;
  const vector<CalibData::ValSig> *psnb;

  CalibData::CalAsym *asym = (CalibData::CalAsym*)(rngBase);

  if (!(lrg = asym->getBig())) {
    msglog << MSG::ERROR << "Unable to retrieve calib data for " << m_calibPath << endreq;
    return false;
  }
  if (!(sm = asym->getSmall())) {
    msglog << MSG::ERROR << "Unable to retrieve calib data for " << m_calibPath << endreq;
    return false;
  }
  if (!(nspb = asym->getNSmallPBig())) {
    msglog << MSG::ERROR << "Unable to retrieve calib data for " << m_calibPath << endreq;
    return false;
  }
  if (!(psnb = asym->getPSmallNBig())) {
    msglog << MSG::ERROR << "Unable to retrieve calib data for " << m_calibPath << endreq;
    return false;
  }

  // get Xpos values.
  unsigned XposSize= m_calibBase->getXpos()->getVals()->size();
  if (XposSize != lrg->size() ||
      XposSize != sm->size() ||
      XposSize != nspb->size() ||
      XposSize != psnb->size()) {
    msglog << MSG::ERROR << "Invalid # of values for " << m_calibPath << endreq;
    return false;
  }
  return true;
}

/// retrieve Asymmetry calibration information for one xtal
StatusCode AsymMgr::getAsym(const CalXtalId &xtalId,
                            const vector<CalibData::ValSig> *&lrg,
                            const vector<CalibData::ValSig> *&sm,
                            const vector<CalibData::ValSig> *&NSmPLrg,
                            const vector<CalibData::ValSig> *&PSmNLrg,
                            const vector<float>  *&xVals) {
  if (!checkXtalId(xtalId)) return StatusCode::FAILURE;

  CalibData::CalAsym *asym = getRangeBase(xtalId);
  if (!asym) return StatusCode::FAILURE;

  // get main data arrays
  lrg = asym->getBig();
  sm = asym->getSmall();
  NSmPLrg = asym->getNSmallPBig();
  PSmNLrg = asym->getPSmallNBig();
  CalibData::Xpos *tmpXpos = m_calibBase->getXpos();
  xVals = tmpXpos->getVals();

  return StatusCode::SUCCESS;
}

StatusCode AsymMgr::genSplines() {
  StatusCode sc;

  // vector<double> arrays for input into genSpline
  vector<double> dblLrg;
  vector<double> dblSm;
  vector<double> dblNSPB;
  vector<double> dblPSNB;
  vector<double> dblXpos;

  for (XtalIdx xtalIdx; xtalIdx.isValid(); xtalIdx++) {
    // needed params for getAsym
    const vector<CalibData::ValSig> *lrg;
    const vector<CalibData::ValSig> *sm;
    const vector<CalibData::ValSig> *nspb;
    const vector<CalibData::ValSig> *psnb;
    const vector<float> *Xpos;

    sc = getAsym(xtalIdx.getCalXtalId(), 
                 lrg, sm, 
                 nspb, psnb, 
                 Xpos);
    if (sc.isFailure()) return sc;

    int n = Xpos->size();
    
    //-- add one point on either end for --//
    //-- linearly extrapolated end-points --//
    dblLrg.resize(n+2);
    dblSm.resize(n+2);
    dblPSNB.resize(n+2);
    dblNSPB.resize(n+2);
    dblXpos.resize(n+2);

    for (int i = 0; i < n; i++) {
      dblLrg[i+1]   = (*lrg)[i].getVal();
      dblSm[i+1] = (*sm)[i].getVal();
      dblNSPB[i+1]  = (*nspb)[i].getVal();
      dblPSNB[i+1]  = (*psnb)[i].getVal();
      dblXpos[i+1]  = (*Xpos)[i];
    }

    //-- LINEAR EXTRAPOLATION
    dblLrg[0]  = 2*dblLrg[1]  - dblLrg[2];
    dblSm[0]   = 2*dblSm[1]   - dblSm[2];
    dblNSPB[0] = 2*dblNSPB[1] - dblNSPB[2];
    dblPSNB[0] = 2*dblPSNB[1] - dblPSNB[2];
    dblXpos[0] = 2*dblXpos[1] - dblXpos[2];

    dblLrg[n+1]  = 2*dblLrg[n]  - dblLrg[n-1];
    dblSm[n+1]   = 2*dblSm[n]   - dblSm[n-1];
    dblNSPB[n+1] = 2*dblNSPB[n] - dblNSPB[n-1];
    dblPSNB[n+1] = 2*dblPSNB[n] - dblPSNB[n-1];
    dblXpos[n+1] = 2*dblXpos[n] - dblPSNB[n-1];

    genSpline(LRG_SPLINE,   xtalIdx, "lrg",   dblXpos, dblLrg);
    genSpline(SM_SPLINE,    xtalIdx, "sm",    dblXpos, dblSm);
    genSpline(NSPB_SPLINE,  xtalIdx, "nspb",  dblXpos, dblNSPB);
    genSpline(PSNB_SPLINE,  xtalIdx, "psnb",  dblXpos, dblPSNB);

    genSpline(INV_LRG_SPLINE,   xtalIdx, "invLrg",   dblLrg,   dblXpos);
    genSpline(INV_SM_SPLINE,    xtalIdx, "invSm",    dblSm,    dblXpos);
    genSpline(INV_NSPB_SPLINE,  xtalIdx, "invNspb",  dblNSPB,  dblXpos);
    genSpline(INV_PSNB_SPLINE,  xtalIdx, "invPsnb",  dblPSNB,  dblXpos);
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
