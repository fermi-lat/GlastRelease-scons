#ifndef AsymMgr_H
#define AsymMgr_H
// $Header$
// LOCAL
#include "CalibItemMgr.h"


// GLAST
#include "CalUtil/CalDefs.h"
#include "CalibData/Cal/Xpos.h"
#include "CalibData/Cal/CalAsym.h"

// EXTLIB

// STD
#include <memory>

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
  CalibData::CalAsym *getAsym(CalUtil::XtalIdx xtalIdx);

  CalibData::Xpos *getXpos();

  StatusCode evalAsym(CalUtil::XtalIdx xtalIdx, CalUtil::AsymType asymType,
                      float Xpos, float &asym){
    return evalSpline(splineIdx(asymType, POS2ASYM), xtalIdx, Xpos, asym);
  }

  StatusCode evalPos(CalUtil::XtalIdx xtalIdx, CalUtil::AsymType asymType,
                     float asym, float &Xpos) {
    return evalSpline(splineIdx(asymType, ASYM2POS), xtalIdx, asym, Xpos);
  }

  StatusCode getAsymCtr(CalUtil::XtalIdx xtalIdx, 
                        CalUtil::AsymType asymType, 
                        float &asymCtr);


 private:
    typedef enum {
      POS2ASYM,
      ASYM2POS,
      N_ASYM_DIR
    } ASYM_DIR;

    static const unsigned short N_SPLINE_TYPES = N_ASYM_DIR*CalUtil::AsymType::N_VALS;

    unsigned short splineIdx(CalUtil::AsymType asymType, ASYM_DIR asymDir) {
      return asymType.val()*N_ASYM_DIR + asymDir;
    }


  StatusCode genLocalStore();

  StatusCode loadIdealVals();

  /// store ideal asymmetry splines (same for every xtal)
  std::auto_ptr<CalibData::CalAsym> m_idealAsym;

  /// store position values for ideal asym splines.
  std::auto_ptr<CalibData::Xpos> m_idealXpos;

  /// store precalcuated asymmetry for deposits @ ctr of xtal.
  /// measures 'overall' optical asymmetry
  CalUtil::CalArray<CalUtil::AsymType, CalUtil::CalArray<CalUtil::XtalIdx, float> > m_asymCtr;
  
  /// validate data ptr from TDS
  bool validateRangeBase(CalibData::CalAsym *asym);

  /// geometry constant
  float m_csiLength;

};
#endif

