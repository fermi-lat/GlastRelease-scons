// $Header$
#include "CalibData/Cal/LightAsym.h"

namespace CalibData {

  LightAsym::LightAsym(const std::vector<float>* values, float err) :
      m_values(0), m_error(err) {
    if (values != 0) {
      m_values = new std::vector<float>(*values);
    }
  }

  void LightAsym::update(RangeBase* other)  {
    LightAsym* otherLightAsym = dynamic_cast<LightAsym* > (other);
    
    // check that otherLightAsym isn't 0 (dynamic cast worked)
    m_values = new std::vector<float>(*otherLightAsym->m_values);
    m_error = otherLightAsym->m_error;
  }
}
