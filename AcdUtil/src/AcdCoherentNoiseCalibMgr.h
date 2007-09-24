#ifndef CoherentNoiseMgr_H
#define CoherentNoiseMgr_H
// $Header$
// LOCAL
#include "AcdCalibMgr.h"

// GLAST
#include "CalibData/Acd/AcdCoherentNoise.h"
#include "CalibData/CalibModel.h"
//#include "CalibSvc/ICalibPathSvc.h"


// EXTLIB
// STD


class AcdCalibSvc;

/** @class AcdCoherentNoiseCalibMgr
    @author Eric Charles (from Zachary Fewtrell's CalCalib stuff)
    
    \brief Manage GLAST Acd pedestal calibration data
*/

class AcdCoherentNoiseCalibMgr : public AcdCalibMgr {
public:
  AcdCoherentNoiseCalibMgr() : 
    //      AcdCalibMgr(ICalibPathSvc::Calib_ACD_CoherentNoise)
    AcdCalibMgr(CalibData::ACD_CoherentNoise)
  {};
  
  /// get pedestal vals for given channel
  StatusCode getCoherentNoise(idents::AcdId id, unsigned pmt,
		      CalibData::AcdCoherentNoise*& calib) {

    static CalibData::AcdCoherentNoise nullCoherentNoise(0.,0.,0.,0.,0); // Amplitude is 0, no oscillation
    if ( m_ideal ) {
      calib = &nullCoherentNoise;
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
    calib = static_cast<CalibData::AcdCoherentNoise*>(base);

    // make sure that it really exists
    if (!calib) return StatusCode::FAILURE;
    
    return StatusCode::SUCCESS;
  }
  
private:

};

#endif
