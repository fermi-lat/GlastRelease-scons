#ifndef VetoMgr_H
#define VetoMgr_H
// $Header$
// LOCAL
#include "AcdCalibMgr.h"

// GLAST
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/CalibModel.h"

// EXTLIB
// STD


class AcdCalibSvc;

/** @class AcdVetoCalibMgr
    @author Eric Charles (from Zachary Fewtrell's CalCalib stuff)
    
    \brief Manage GLAST Acd pedestal calibration data
*/

class AcdVetoCalibMgr : public AcdCalibMgr {
public:
  AcdVetoCalibMgr() : 
    AcdCalibMgr(CalibData::ACD_ThreshVeto)
  {};
  
  /// get pedestal vals for given channel
  StatusCode getVeto(idents::AcdId id, unsigned pmt,
		    CalibData::AcdVeto*& veto) {

    static CalibData::AcdVeto nullVeto(200.,0.,0); // veto fires at 200 counts PHA
    if ( m_ideal ) {
      veto = &nullVeto;
      return StatusCode::SUCCESS;
    }
  
    // make sure we have valid calib data for this event.
    StatusCode sc;
    sc = updateCalib();
    if (sc.isFailure()) {
      // null and return failure code
      veto = 0;
      return sc;
    }
 
    // grab the data we need
    CalibData::RangeBase* base = m_calibBase->getPmt(id,pmt);    
    veto = static_cast<CalibData::AcdVeto*>(base);

    // make sure that it really exists
    if (!veto) return StatusCode::FAILURE;
    
    return StatusCode::SUCCESS;
  }
  
private:
  
};

#endif
