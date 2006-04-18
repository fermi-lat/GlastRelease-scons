#ifndef PedMgr_H
#define PedMgr_H
// $Header$
// LOCAL
#include "AcdCalibMgr.h"

// GLAST
#include "CalibData/Acd/AcdPed.h"
#include "CalibData/CalibModel.h"

// EXTLIB
// STD


class AcdCalibSvc;

/** @class AcdPedCalibMgr
    @author Eric Charles (from Zachary Fewtrell's CalCalib stuff)
    
    \brief Manage GLAST Acd pedestal calibration data
*/

class AcdPedCalibMgr : public AcdCalibMgr {
public:
  AcdPedCalibMgr() : 
    AcdCalibMgr(CalibData::ACD_Ped)
  {};
  
  /// get pedestal vals for given channel
  StatusCode getPed(idents::AcdId id, unsigned pmt,
		    CalibData::AcdPed*& pedestal) {

    static CalibData::AcdPed nullPed(0.,0.,0); // all pedestals are null
    if ( m_ideal ) {
      pedestal = &nullPed;
      return StatusCode::SUCCESS;
    }
  
    // make sure we have valid calib data for this event.
    StatusCode sc;
    sc = updateCalib();
    if (sc.isFailure()) {
      // null and return failure code
      pedestal = 0;
      return sc;
    }
 
    // grab the data we need
    CalibData::RangeBase* base = m_calibBase->getPmt(id,pmt);    
    pedestal = static_cast<CalibData::AcdPed*>(base);

    // make sure that it really exists
    if (!pedestal) return StatusCode::FAILURE;
    
    return StatusCode::SUCCESS;
  }
  
private:
  
};

#endif
