#ifndef AsymMgr_H
#define AsymMgr_H 1

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalibData/Cal/CalAsym.h"

// EXTLIB

// STD

using namespace CalDefs;
using namespace idents;

class AsymMgr : public CalibItemMgr {
 public:
  AsymMgr();

  /// retrieve Asymmetry calibration information for one xtal
  StatusCode getAsym(const CalXtalId &xtalId,
                     const vector<CalibData::ValSig> *&lrg,
                     const vector<CalibData::ValSig> *&sm,
                     const vector<CalibData::ValSig> *&NSmPLrg,
                     const vector<CalibData::ValSig> *&PSmNLrg,
                     const vector<float>  *&xVals);

  StatusCode evalAsymLrg(const CalXtalId &xtalId, double Xpos, double &lrg) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(LRG_SPLINE, XtalIdx(xtalId), Xpos, lrg);
  }

  StatusCode evalPosLrg(const CalXtalId &xtalId, double lrg, double &Xpos) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(INV_LRG_SPLINE, XtalIdx(xtalId), lrg, Xpos);
  }

  StatusCode evalAsymSm(const CalXtalId &xtalId, double Xpos, double &sm) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(SM_SPLINE, XtalIdx(xtalId), Xpos, sm);
  }

  StatusCode evalPosSm(const CalXtalId &xtalId, double sm, double &Xpos) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(INV_SM_SPLINE, XtalIdx(xtalId), sm, Xpos);
  }

  StatusCode evalAsymNSPB(const CalXtalId &xtalId, double Xpos, double &nspb) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(NSPB_SPLINE, XtalIdx(xtalId), Xpos, nspb);
  }

  StatusCode evalPosNSPB(const CalXtalId &xtalId, double nspb, double &Xpos) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(INV_NSPB_SPLINE, XtalIdx(xtalId), nspb, Xpos);
  }

  StatusCode evalAsymPSNB(const CalXtalId &xtalId, double Xpos, double &psnb) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(PSNB_SPLINE, XtalIdx(xtalId), Xpos, psnb);
  }

  StatusCode evalPosPSNB(const CalXtalId &xtalId, double psnb, double &Xpos) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(INV_PSNB_SPLINE, XtalIdx(xtalId), psnb, Xpos);
  }

 private:

  typedef enum {
    LRG_SPLINE,
    INV_LRG_SPLINE,
    SM_SPLINE,
    INV_SM_SPLINE,
    NSPB_SPLINE,
    INV_NSPB_SPLINE,
    PSNB_SPLINE,
    INV_PSNB_SPLINE,
    N_SPLINE_TYPES
  } SPLINE_TYPE;

  CalibData::CalAsym *getRangeBase(const CalXtalId &xtalId) {
    StatusCode sc = updateCache();
    if (sc.isFailure()) return NULL;
    return (CalibData::CalAsym*)
      (m_rngBases[XtalIdx(xtalId)]);
  }

  bool validateRangeBase(const CalXtalId &xtalId, CalibData::RangeBase *rngBase);

  StatusCode fillRangeBases();

  StatusCode genSplines();

  bool checkXtalId(const CalXtalId &xtalId) {return true;}
};
#endif
