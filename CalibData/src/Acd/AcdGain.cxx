// $Header$
#include "CalibData/Acd/AcdGain.h"

namespace CalibData {

  void AcdGain::update(RangeBase* other)  {
    AcdGain* otherGain = dynamic_cast<AcdGain* > (other);
    
    // check that otherGain isn't 0 (dynamic cast worked)
    m_peak = otherGain->m_peak;
    m_width = otherGain->m_width;
    m_status = otherGain->m_status;
  }
}
