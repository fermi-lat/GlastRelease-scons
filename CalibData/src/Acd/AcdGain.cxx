// $Header$
#include "CalibData/Acd/AcdGain.h"

namespace CalibData {

  void AcdGain::update(RangeBase* other)  {
    AcdGain* otherGain = dynamic_cast<AcdGain* > (other);
    
    // check that otherGain isn't 0 (dynamic cast worked)
    m_gain = otherGain->m_gain;
    m_sig = otherGain->m_sig;
  }
}
