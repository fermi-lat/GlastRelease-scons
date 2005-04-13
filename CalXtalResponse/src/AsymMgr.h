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

using CalibData::ValSig;

class CalCalibSvc;

class AsymMgr : public CalibItemMgr {
 public:
  AsymMgr();

  /// get Asymmetry calibration information for one xtal
  StatusCode getAsym(const CalXtalId &xtalId,
                     const vector<ValSig> *&asymLrg,
                     const vector<ValSig> *&asymSm,
                     const vector<ValSig> *&asymNSPB,
                     const vector<ValSig> *&asymPSNB,
                     const vector<float>  *&xVals);

  StatusCode evalAsymLrg(const CalXtalId &xtalId, 
                         double Xpos, double &asymLrg){
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(ASYMLRG_SPLINE, XtalIdx(xtalId), Xpos, asymLrg);
  }

  StatusCode evalPosLrg(const CalXtalId &xtalId, 
                        double asymLrg, double &Xpos) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(INV_ASYMLRG_SPLINE, XtalIdx(xtalId), asymLrg, Xpos);
  }

  StatusCode evalAsymSm(const CalXtalId &xtalId, 
                        double Xpos, double &asymSm) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(ASYMSM_SPLINE, XtalIdx(xtalId), Xpos, asymSm);
  }

  StatusCode evalPosSm(const CalXtalId &xtalId, 
                       double asymSm, double &Xpos) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(INV_ASYMSM_SPLINE, XtalIdx(xtalId), asymSm, Xpos);
  }

  StatusCode evalAsymNSPB(const CalXtalId &xtalId, 
                          double Xpos, double &asymNSPB) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(ASYMNSPB_SPLINE, XtalIdx(xtalId), Xpos, asymNSPB);
  }

  StatusCode evalPosNSPB(const CalXtalId &xtalId, 
                         double asymNSPB, double &Xpos) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(INV_ASYMNSPB_SPLINE, XtalIdx(xtalId), asymNSPB, Xpos);
  }

  StatusCode evalAsymPSNB(const CalXtalId &xtalId, 
                          double Xpos, double &asymPSNB) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(ASYMPSNB_SPLINE, XtalIdx(xtalId), Xpos, asymPSNB);
  }

  StatusCode evalPosPSNB(const CalXtalId &xtalId, 
                         double asymPSNB, double &Xpos) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(INV_ASYMPSNB_SPLINE, XtalIdx(xtalId), asymPSNB, Xpos);
  }

 private:

  /// List of each type of spline for asym calib_type
  typedef enum {
    ASYMLRG_SPLINE,
    INV_ASYMLRG_SPLINE,
    ASYMSM_SPLINE,
    INV_ASYMSM_SPLINE,
    ASYMNSPB_SPLINE,
    INV_ASYMNSPB_SPLINE,
    ASYMPSNB_SPLINE,
    INV_ASYMPSNB_SPLINE,
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

  bool checkXtalId(const CalXtalId&) {return true;}

  StatusCode loadIdealVals();
  /// Store ideal (fake) vals for large diode asym (used when db is down)
  vector<ValSig> m_idealAsymLrg;
  /// Store ideal (fake) vals for small diode asym (used when db is down)
  vector<ValSig> m_idealAsymSm;
  /// Store ideal (fake) vals for NegSmall diode PosBig diode asym (used when db is down)
  vector<ValSig> m_idealAsymNSPB;
  /// Store ideal (fake) vals for PosSmall diode NegBig diode asym (used when db is down)
  vector<ValSig> m_idealAsymPSNB;
  
  vector<float> m_idealXVals;

};
#endif
