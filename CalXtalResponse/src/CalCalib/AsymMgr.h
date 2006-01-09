#ifndef AsymMgr_H
#define AsymMgr_H
// $Header$
// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalibData/Cal/CalAsym.h"

// EXTLIB

// STD

using namespace CalUtil;
using namespace idents;

using CalibData::ValSig;

class CalCalibSvc;

/** @class AsymMgr
    @author Zachary Fewtrell
    
    \brief Manage GLAST Cal asymmetry calibration data
*/

class AsymMgr : public CalibItemMgr {
 public:
  AsymMgr();

  /// get Asymmetry calibration information for one xtal
  StatusCode getAsym(XtalIdx xtalIdx,
                     const vector<ValSig> *&asymLrg,
                     const vector<ValSig> *&asymSm,
                     const vector<ValSig> *&asymNSPB,
                     const vector<ValSig> *&asymPSNB,
                     const vector<float>  *&xVals);

  StatusCode evalAsymLrg(XtalIdx xtalIdx, 
                         float Xpos, float &asymLrg){
    return evalSpline(ASYMLRG_SPLINE, xtalIdx, Xpos, asymLrg);
  }

  StatusCode evalPosLrg(XtalIdx xtalIdx, 
                        float asymLrg, float &Xpos) {
    return evalSpline(INV_ASYMLRG_SPLINE, xtalIdx, asymLrg, Xpos);
  }

  StatusCode evalAsymSm(XtalIdx xtalIdx, 
                        float Xpos, float &asymSm) {
    return evalSpline(ASYMSM_SPLINE, xtalIdx, Xpos, asymSm);
  }

  StatusCode evalPosSm(XtalIdx xtalIdx, 
                       float asymSm, float &Xpos) {
    return evalSpline(INV_ASYMSM_SPLINE, xtalIdx, asymSm, Xpos);
  }

  StatusCode evalAsymNSPB(XtalIdx xtalIdx, 
                          float Xpos, float &asymNSPB) {
    return evalSpline(ASYMNSPB_SPLINE, xtalIdx, Xpos, asymNSPB);
  }

  StatusCode evalPosNSPB(XtalIdx xtalIdx, 
                         float asymNSPB, float &Xpos) {
    return evalSpline(INV_ASYMNSPB_SPLINE, xtalIdx, asymNSPB, Xpos);
  }

  StatusCode evalAsymPSNB(XtalIdx xtalIdx, 
                          float Xpos, float &asymPSNB) {
    return evalSpline(ASYMPSNB_SPLINE, xtalIdx, Xpos, asymPSNB);
  }

  StatusCode evalPosPSNB(XtalIdx xtalIdx, 
                         float asymPSNB, float &Xpos) {
    return evalSpline(INV_ASYMPSNB_SPLINE, xtalIdx, asymPSNB, Xpos);
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

  StatusCode genLocalStore();

  StatusCode loadIdealVals();
  /// Store ideal (fake) vals for large diode asym (used when db is down)
  vector<ValSig> m_idealAsymLrg;
  /// Store ideal (fake) vals for small diode asym (used when db is down)
  vector<ValSig> m_idealAsymSm;
  /// Store ideal (fake) vals for NegSmall diode PosBig diode asym 
  /// (used when db is down)

  vector<ValSig> m_idealAsymNSPB;
  /// Store ideal (fake) vals for PosSmall diode NegBig diode asym 
  /// (used when db is down)
  vector<ValSig> m_idealAsymPSNB;
  
  vector<float> m_idealXVals;

  /// validate data ptr from TDS
  bool validateRangeBase(CalibData::CalAsym *asym);

};
#endif
