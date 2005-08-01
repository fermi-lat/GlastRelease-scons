// $Header$
#include "CalibData/Cal/IntNonlin.h"

namespace CalibData {

  IntNonlin::IntNonlin(const std::vector<float>* values, float err, 
                       const std::vector<float>* sdacs) :
      m_values(0), m_sdacs(0), m_error(err) {
    if (values != 0) {
      m_values = new std::vector<float>(*values);
    }
    if (sdacs != 0) {
      m_sdacs = new std::vector<float>(*sdacs);
    }
  }

  void IntNonlin::update(RangeBase* other)  {
    IntNonlin* otherIntNonlin = dynamic_cast<IntNonlin* > (other);
    
    // check that otherIntNonlin isn't 0 (dynamic cast worked)
    m_values = new std::vector<float>(*otherIntNonlin->m_values);
    if (otherIntNonlin->m_sdacs) {
      m_sdacs = new std::vector<float>(*otherIntNonlin->m_sdacs);
    }
    m_error = otherIntNonlin->m_error;
  }
}
