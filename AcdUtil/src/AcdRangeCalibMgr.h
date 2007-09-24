#ifndef RangeMgr_H
#define RangeMgr_H
// $Header$
// LOCAL
#include "AcdCalibMgr.h"

// GLAST
#include "CalibData/Acd/AcdRange.h"
#include "CalibData/CalibModel.h"
//#include "CalibSvc/ICalibPathSvc.h"


// EXTLIB
// STD


class AcdCalibSvc;

/** @class AcdRangeCalibMgr
    @author Eric Charles (from Zachary Fewtrell's CalCalib stuff)
    
    \brief Manage GLAST Acd pedestal calibration data
*/

class AcdRangeCalibMgr : public AcdCalibMgr {
public:
  AcdRangeCalibMgr() : 
    // AcdCalibMgr(ICalibPathSvc::Calib_ACD_Range)
    AcdCalibMgr(CalibData::ACD_Range)
  {};
  
  /// get pedestal vals for given channel
  StatusCode getRange(idents::AcdId id, unsigned pmt,
		      CalibData::AcdRange*& range) {

    static CalibData::AcdRange nullRange(4095.,0.,0); // Switch occurs at 4095 in low range = 0 in High Range
    if ( m_ideal ) {
      range = &nullRange;
      return StatusCode::SUCCESS;
    }
  
    // make sure we have valid calib data for this event.
    StatusCode sc;
    sc = updateCalib();
    if (sc.isFailure()) {
      // null and return failure code
      range = 0;
      return sc;
    }
 
    // grab the data we need
    CalibData::RangeBase* base = m_calibBase->getPmt(id,pmt);    
    range = static_cast<CalibData::AcdRange*>(base);

    // make sure that it really exists
    if (!range) return StatusCode::FAILURE;
    
    return StatusCode::SUCCESS;
  }
  
private:

};

#endif
