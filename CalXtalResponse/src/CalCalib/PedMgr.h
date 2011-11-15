#ifndef PedMgr_H
#define PedMgr_H
// $Header$
/** @file 
    @author Z.Fewtrell
*/

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"
#include "CalUtil/CalArray.h"
#include "CalibData/Cal/Ped.h"

// EXTLIB
// STD

class CalCalibSvc;

/** @class PedMgr
    @author Z.Fewtrell
    
    \brief Manage GLAST Cal pedestal calibration data
*/

class PedMgr : public CalibItemMgr {
 public:
  PedMgr(CalCalibShared &ccsShared) : 
    CalibItemMgr(ICalibPathSvc::Calib_CAL_Ped, 
                 ccsShared,
                 CalUtil::RngIdx::N_VALS
                 )
    {};

  /// get pedestal vals for given xtal/face/rng
  const CalibData::Ped *getPed(const CalUtil::RngIdx rngIdx);
 private:
  /// load ideal calibration data from local store (bypass calib db)
  StatusCode loadIdealVals();

  /// validate calibration data for single channel
  bool validateRangeBase(CalibData::Ped const*const ped);
  
  /// ideal calibration data (same for all xtals)
  CalUtil::CalArray<CalUtil::RngNum, CalibData::Ped> m_idealPeds;

  /// populate all internal calibration data fields
  StatusCode genLocalStore();
};

#endif
