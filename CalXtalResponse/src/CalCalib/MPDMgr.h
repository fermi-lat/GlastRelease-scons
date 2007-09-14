#ifndef MPDMgr_H
#define MPDMgr_H
// $Header$
// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/CalMevPerDac.h"

// EXTLIB

// STD
#include <memory>

class CalCalibSvc;

/** \brief Manage MevPerDAC calibration type
    \author Zach Fewtrell
    */
class MPDMgr : public CalibItemMgr {
 public:
  MPDMgr(CalCalibShared &ccsShared) :
    CalibItemMgr(ICalibPathSvc::Calib_CAL_MevPerDac,
                 ccsShared) {}

  const CalibData::CalMevPerDac *getMPD(const CalUtil::XtalIdx xtalIdx);

 private:
  StatusCode loadIdealVals();
  
  StatusCode genLocalStore();

  static bool validateRangeBase(const CalibData::CalMevPerDac *mpd);

  /// ideal data is not available till after construction & I will be responsible
  /// for memory, hence use of auto_ptr<>
  std::auto_ptr<CalibData::CalMevPerDac> m_idealMPD;
};

#endif
