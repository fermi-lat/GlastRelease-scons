// $Header$
#include "CalibData/Cal/Gain.h"

namespace CalibData {

  void Gain::update(RangeBase* other)  {
    Gain* otherGain = dynamic_cast<Gain* > (other);
    
    // check that otherGain isn't 0 (dynamic cast worked)
    m_gain = otherGain->m_gain;
    m_sig = otherGain->m_sig;
  }
}
