#ifndef PedMgr_H
#define PedMgr_H
// $Header$
// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"
#include "CalibData/Cal/Ped.h"

// EXTLIB
// STD

class CalCalibSvc;

/** @class PedMgr
    @author Zachary Fewtrell
    
    \brief Manage GLAST Cal pedestal calibration data
*/

class PedMgr : public CalibItemMgr {
 public:
  PedMgr(CalCalibShared &ccsShared) : 
    CalibItemMgr(ICalibPathSvc::Calib_CAL_Ped, 
                 ccsShared)
    {};

  /// get pedestal vals for given xtal/face/rng
  const CalibData::Ped *getPed(CalUtil::RngIdx rngIdx);
 private:
  StatusCode loadIdealVals();

  bool validateRangeBase(CalibData::Ped *ped);
  
  CalUtil::CalArray<CalUtil::RngNum, CalibData::Ped> m_idealPeds;

  StatusCode genLocalStore();
};

#endif
