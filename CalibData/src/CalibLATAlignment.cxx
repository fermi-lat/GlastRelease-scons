// $Header$

/** @class CalibLATAlignment
 *    Implementation of a calibration class for the LAT-to-spacecraft alignment
 */

#include "CalibData/CalibLATAlignment.h"
#include "GaudiKernel/MsgStream.h"

namespace CalibData {
  
  CalibLATAlignment::CalibLATAlignment(double px, double py, double pz,
                                        const ITime& since, const ITime& till,int serNo):
         CalibBase(since, till, serNo), m_roll(px), m_pitch(py), m_yaw(pz) {
    //    m_me = this;
 
  }  

  CalibLATAlignment::CalibLATAlignment(const CalibLATAlignment& other) : 
        IValidity(), CalibBase(other), m_roll(other.m_roll), m_pitch(other.m_pitch), m_yaw(other.m_yaw){
  }

  StatusCode CalibLATAlignment::update(CalibBase& other, MsgStream* log) {
    // The following dynamic_cast has got to work
    CalibLATAlignment& other1 = dynamic_cast<CalibLATAlignment& >(other);

    CalibBase::update(other1, log);
    m_roll = other1.m_roll;
    m_pitch = other1.m_pitch;
    m_yaw = other1.m_yaw;

    return StatusCode::SUCCESS;
  }


}
