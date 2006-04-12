#ifndef AcdGainCalibMgr_H
#define AcdGainCalibMgr_H
// $Header$
// LOCAL
#include "AcdCalibMgr.h"

// GLAST
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/CalibModel.h"

// EXTLIB
// STD

class AcdCalibSvc;

/** @class AcdGainCalibMgr
    @author Eric Charles (from Zachary Fewtrell's CalCalib stuff)
    
    \brief Manage GLAST Cal pedestal calibration data
*/

class AcdGainCalibMgr : public AcdCalibMgr {
public:
  AcdGainCalibMgr() : 
    AcdCalibMgr(CalibData::ACD_ElecGain)
  {};
  
  /// get mip peak vals for given channel
  StatusCode getMipPeak(idents::AcdId id, unsigned pmt,
			CalibData::AcdGain*& mipPeak) {

    // make sure we have valid calib data for this event.
    StatusCode sc;
    sc = updateCalib();
    if (sc.isFailure()){
      // null and return failure code
      mipPeak = 0;
      return sc;
    }
    
    // grab the data we need
    CalibData::RangeBase* base = m_calibBase->getPmt(id,pmt);    
    mipPeak = static_cast<CalibData::AcdGain*>(base);

    // make sure that it really exists
    if (!mipPeak) return StatusCode::FAILURE;
    
    return StatusCode::SUCCESS;
  }
  
private:
  
};

#endif
