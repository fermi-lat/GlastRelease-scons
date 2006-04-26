#ifndef AsymMgr_H
#define AsymMgr_H
// $Header$
// LOCAL
#include "CalibItemMgr.h"
#include "CalXtalResponse/CalCalibDefs.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalibData/Cal/CalAsym.h"
#include "CalibData/Cal/Xpos.h"

// EXTLIB

// STD

using namespace CalUtil;
using namespace idents;
using namespace CalXtalResponse;
using namespace CalibData;

class CalCalibSvc;

/** @class AsymMgr
    @author Zachary Fewtrell
    
    \brief Manage GLAST Cal asymmetry calibration data
*/

class AsymMgr : public CalibItemMgr {
 public:
  AsymMgr(CalCalibShared &ccsShared);

  StatusCode initialize(const string &flavor);

  /// get Asymmetry calibration information for one xtal
  CalAsym *getAsym(XtalIdx xtalIdx);

  Xpos *getXpos();

  StatusCode evalAsym(XtalIdx xtalIdx, AsymType asymType,
                      float Xpos, float &asym){
    return evalSpline(splineIdx(asymType, POS2ASYM), xtalIdx, Xpos, asym);
  }

  StatusCode evalPos(XtalIdx xtalIdx, AsymType asymType,
                     float asym, float &Xpos) {
    return evalSpline(splineIdx(asymType, ASYM2POS), xtalIdx, asym, Xpos);
  }

  StatusCode getAsymCtr(XtalIdx xtalIdx, 
                        AsymType asymType, 
                        float &asymCtr);


 private:
    typedef enum {
      POS2ASYM,
      ASYM2POS,
      N_ASYM_DIR
    } ASYM_DIR;

    static const unsigned short N_SPLINE_TYPES = N_ASYM_DIR*AsymType::N_VALS;

    unsigned short splineIdx(AsymType asymType, ASYM_DIR asymDir) {
      return asymType.val()*N_ASYM_DIR + asymDir;
    }


  StatusCode genLocalStore();

  StatusCode loadIdealVals();

  /// store ideal asymmetry splines (same for every xtal)
  auto_ptr<CalAsym> m_idealAsym;

  /// store position values for ideal asym splines.
  auto_ptr<Xpos> m_idealXpos;

  /// store precalcuated asymmetry for deposits @ ctr of xtal.
  /// measures 'overall' optical asymmetry
  CalArray<AsymType, CalArray<XtalIdx, float> > m_asymCtr;
  
  /// validate data ptr from TDS
  bool validateRangeBase(CalAsym *asym);

  /// geometry constant
  float m_csiLength;

};
#endif

