// $Header$

/** @class CalibTest1
 *    Implementation of near-simplest-possible calibration TCDS class
 */

#include "CalibData/CalibTest1.h"

namespace CalibData {
  void CalibTest1::update(CalibTest1& other) {
    CalibBase::update(other);
    m_name = other.m_name;
    m_value = other.m_value;
  }

  CalibTest1::CalibTest1(const CalibTest1& other) : CalibBase(other),
                                              m_name(other.m_name), 
                                              m_value(other.m_value){}
    
  
  
}
