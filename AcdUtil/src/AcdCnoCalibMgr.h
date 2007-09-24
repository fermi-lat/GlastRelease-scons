#ifndef CnoMgr_H
#define CnoMgr_H
// $Header$
// LOCAL
#include "AcdCalibMgr.h"

// GLAST
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/CalibModel.h"

// EXTLIB
// STD


class AcdCalibSvc;

/** @class AcdCnoCalibMgr
    @author Eric Charles (from Zachary Fewtrell's CalCalib stuff)
    
    \brief Manage GLAST Acd pedestal calibration data
*/

class AcdCnoCalibMgr : public AcdCalibMgr {
public:
  AcdCnoCalibMgr() : 
    AcdCalibMgr(CalibData::ACD_ThreshHigh)
  {};
  
  /// get pedestal vals for given channel
  StatusCode getCno(idents::AcdId id, unsigned pmt,
		    CalibData::AcdCno*& cno) {

    static CalibData::AcdCno nullCno(0.,0.,0); // all pedestals are null
    if ( m_ideal ) {
      cno = &nullCno;
      return StatusCode::SUCCESS;
    }
  
    // make sure we have valid calib data for this event.
    StatusCode sc;
    sc = updateCalib();
    if (sc.isFailure()) {
      // null and return failure code
      cno = 0;
      return sc;
    }
 
    // grab the data we need
    CalibData::RangeBase* base = m_calibBase->getPmt(id,pmt);    
    cno = static_cast<CalibData::AcdCno*>(base);

    // make sure that it really exists
    if (!cno) return StatusCode::FAILURE;
    
    return StatusCode::SUCCESS;
  }
  
private:
  
};

#endif
