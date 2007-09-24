#ifndef HighRangeMgr_H
#define HighRangeMgr_H
// $Header$
// LOCAL
#include "AcdCalibMgr.h"

// GLAST
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/CalibModel.h"
//#include "CalibSvc/ICalibPathSvc.h"


// EXTLIB
// STD


class AcdCalibSvc;

/** @class AcdHighRangeCalibMgr
    @author Eric Charles (from Zachary Fewtrell's CalCalib stuff)
    
    \brief Manage GLAST Acd pedestal calibration data
*/

class AcdHighRangeCalibMgr : public AcdCalibMgr {
public:
  AcdHighRangeCalibMgr() : 
    //      AcdCalibMgr(ICalibPathSvc::Calib_ACD_HighRange)
    AcdCalibMgr(CalibData::ACD_HighRange)
  {};
  
  /// get pedestal vals for given channel
  StatusCode getHighRange(idents::AcdId id, unsigned pmt,
			  CalibData::AcdHighRange*& calib) {

    static CalibData::AcdHighRange nullHighRange(0.,0.,0); // all pedestals are null
    if ( m_ideal ) {
      calib = &nullHighRange;
      return StatusCode::SUCCESS;
    }
  
    // make sure we have valid calib data for this event.
    StatusCode sc;
    sc = updateCalib();
    if (sc.isFailure()) {
      // null and return failure code
      calib = 0;
      return sc;
    }
 
    // grab the data we need
    CalibData::RangeBase* base = m_calibBase->getPmt(id,pmt);    
    calib = static_cast<CalibData::AcdHighRange*>(base);

    // make sure that it really exists
    if (!calib) return StatusCode::FAILURE;
    
    return StatusCode::SUCCESS;
  }
  
private:

};

#endif
