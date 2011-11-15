#ifndef AsymMgr_H
#define AsymMgr_H
// $Header$
/** @file
    @author Z.Fewtrell
*/
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
    @author Z.Fewtrell
    
    \brief Manage GLAST Cal asymmetry calibration data
*/

class AsymMgr : public CalibItemMgr {
public:
  AsymMgr(CalCalibShared &ccsShared);

  StatusCode initialize(const string &flavor);

  /// get Asymmetry calibration information for one xtal
  const CalibData::CalAsym *getAsym(const CalUtil::XtalIdx xtalIdx);

  /// get positions in mm from crystal center of each asymmetry data point
  const CalibData::Xpos *getXpos();

  /// find light asymmetry value associated with given crystal and
  /// longitudinal position
  StatusCode evalAsym(const CalUtil::XtalIdx xtalIdx, 
                      const CalUtil::AsymType asymType,
                      const float Xpos, 
                      float &asym) {
    return evalSpline(splineIdx(asymType, POS2ASYM), xtalIdx, Xpos, asym);
  }

  /// find crystal longitudinal position (mm from xtal center) associated with given signal asymmetry
  StatusCode evalPos(const CalUtil::XtalIdx xtalIdx, 
                     const CalUtil::AsymType asymType,
                     const float asym, 
                     float &Xpos) {
    return evalSpline(splineIdx(asymType, ASYM2POS), xtalIdx, asym, Xpos);
  }

  /// get signal asymmetry at crystal center.
  StatusCode getAsymCtr(const CalUtil::XtalIdx xtalIdx, 
                        const CalUtil::AsymType asymType, 
                        float &asymCtr);


private:

  /// conversion directions between asymmetry & position.
  typedef enum {
    POS2ASYM,
    ASYM2POS,
    N_ASYM_DIR
  } ASYM_DIR;

  static const unsigned short N_SPLINE_TYPES = N_ASYM_DIR*CalUtil::AsymType::N_VALS;

  /// calculate spline index for given spline type
  unsigned short splineIdx(const CalUtil::AsymType asymType, 
                           const ASYM_DIR asymDir) {
    return asymType.val()*N_ASYM_DIR + asymDir;
  }

  /// populate all private data fields
  StatusCode genLocalStore();

  /// load ideal calibration constants from local data (avoid
  /// calibration database)
  StatusCode loadIdealVals();

  /// store ideal asymmetry splines (same for every xtal)
  std::auto_ptr<CalibData::CalAsym> m_idealAsym;

  /// store position values for ideal asym splines.
  std::auto_ptr<CalibData::Xpos> m_idealXpos;

  /// store precalcuated asymmetry for deposits @ ctr of xtal.
  /// measures 'overall' optical asymmetry
  CalUtil::CalVec<CalUtil::AsymType, CalUtil::CalVec<CalUtil::XtalIdx, float> > m_asymCtr;
  
  /// validate data ptr from TDS
  bool validateRangeBase(CalibData::CalAsym const*const asym);

  /// geometry constant
  float m_csiLength;

};
#endif
