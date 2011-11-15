#ifndef TholdCIMgr_H
#define TholdCIMgr_H
// $Header$
/** @file 
    @author Z.Fewtrell
*/
// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalibData/Cal/CalTholdCI.h"

// EXTLIB

// STD
#include <memory>


class CalCalibSvc;

/** @class TholdCIMgr
    @author Z.Fewtrell
    
    \brief Manage GLAST Cal charge-injection measured threshold data.
*/

class TholdCIMgr : public CalibItemMgr {
 public:
  TholdCIMgr(CalCalibShared &ccsShared) 
    : CalibItemMgr(ICalibPathSvc::Calib_CAL_TholdCI, 
                   ccsShared,
                   CalUtil::FaceIdx::N_VALS
                   )
    {};

  const CalibData::CalTholdCI *getTholdCI(const CalUtil::FaceIdx faceIdx);

 private:
  /// load ideal calibration data from local store (bypass calib db)
  StatusCode loadIdealVals();

  /// populate all internal calibration dat fields
  StatusCode genLocalStore();
  
  /// validate calibration data for single crystal
  bool validateRangeBase(CalibData::CalTholdCI const *const tholdCI);

  /// ideal calibration data (same for all xtals)
  std::auto_ptr<CalibData::CalTholdCI> m_idealTholdCI;
  
};

#endif
