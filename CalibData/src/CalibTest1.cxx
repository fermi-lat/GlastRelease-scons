// $Header$

/** @class CalibTest1
 *    Implementation of near-simplest-possible calibration TCDS class
 */

#include "CalibData/CalibTest1.h"

namespace CalibData {
  CalibTest1::CalibTest1(const std::string& name, int value, 
                         const ITime& since, const ITime& till, 
                         int serNo) :
    CalibBase(since, till, serNo), m_name(name), m_value(value) {
  }


  void CalibTest1::update(CalibTest1& other) {
    CalibBase::update(other);
    m_name = other.m_name;
    m_value = other.m_value;
  }

  CalibTest1::CalibTest1(const CalibTest1& other) : CalibBase(other),
                                                    m_name(other.m_name), 
                                                    m_value(other.m_value) {
  }
    
  std::string CalibTest1::getValueName() const {
    return m_name;
  }
  
}
