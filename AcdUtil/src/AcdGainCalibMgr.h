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

    // four kinds of channels
    static CalibData::AcdGain ribbonGain(56.875,25.4,0);   // ribbons
    static CalibData::AcdGain tileGain(204.75,50.,0);      // most tiles
    static CalibData::AcdGain tile_12mmGain(245.7,50.,0);  // 12mm thick tiles
    static CalibData::AcdGain naGain(-1.,0.,0);            // NA channels
    if ( m_ideal ) {
      if ( id.ribbon() ) {
	mipPeak = &ribbonGain;
      } else if ( id.tile() ) {
	if ( id.face() == 0 && id.row() == 2 ) {
	  mipPeak = &tile_12mmGain;
	} else {
	  mipPeak = &tileGain;
	}
      } else {
	mipPeak = &naGain;
      }
      return StatusCode::SUCCESS;
    }

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
