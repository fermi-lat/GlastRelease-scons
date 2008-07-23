// $Header$

/** @class CalibLATAlignment
 *    Implementation of a calibration class for the LAT-to-spacecraft alignment
 */

#include "CalibData/Nas/CalibLATAlignment.h"
#include "GaudiKernel/MsgStream.h"

namespace CalibData {
  
  CalibLATAlignment::CalibLATAlignment(double px, double py, double pz,
                                       const std::string& units,
                                       const ITime& since, const ITime& till,int serNo):
    CalibBase(since, till, serNo), m_units(units) {
    m_r[0] = px; m_r[1] = py; m_r[2] = pz;
    //    m_me = this;
   }  

  CalibLATAlignment::CalibLATAlignment(const CalibLATAlignment& other) : 
    IValidity(), CalibBase(other), m_units(other.m_units) {
    m_r[0] =other.m_r[0];
    m_r[1] =other.m_r[1];
    m_r[2] =other.m_r[2];
  }

  StatusCode CalibLATAlignment::update(CalibBase& other, MsgStream* log) {
    // The following dynamic_cast has got to work
    CalibLATAlignment& other1 = dynamic_cast<CalibLATAlignment& >(other);

    CalibBase::update(other1, log);
    m_r[0] = other1.m_r[0];
    m_r[1] = other1.m_r[1];
    m_r[2] = other1.m_r[2];
    m_units = other1.m_units;

    return StatusCode::SUCCESS;
  }


}
