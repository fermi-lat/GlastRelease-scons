// $Header$

/** @class CalibTest1
 *    Implementation of near-simplest-possible calibration TCDS class
 */

#include "CalibData/CalibTest1.h"

namespace CalibData {
  CalibTest1::CalibTest1(const std::string& name, int value, 
                         const ITime& since, const ITime& till, 
                         int serNo) :
    CalibBase(since, till, serNo), m_name(name), m_value(value)
  {
    //    m_me = this;
  }


  void CalibTest1::update(CalibBase& other) {
    // The following dynamic_cast has got to work
    CalibTest1& other1 = dynamic_cast<CalibTest1& >(other);

    CalibBase::update(other1);
    m_name = other1.m_name;
    m_value = other1.m_value;
  }

  CalibTest1::CalibTest1(const CalibTest1& other) : CalibBase(other),
                                                    m_name(other.m_name), 
                                                    m_value(other.m_value) {
  }
    
  std::string CalibTest1::getValueName() const {
    return m_name;
  }
  
}
