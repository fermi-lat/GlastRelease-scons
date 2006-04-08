// $Header$
#include "CalibData/Acd/AcdPed.h"

namespace CalibData {

  void AcdPed::update(RangeBase* other)  {
    AcdPed* otherPed = dynamic_cast<AcdPed* > (other);
    
    // check that otherPed isn't 0 (dynamic cast worked)
    m_mean = otherPed->m_mean;
    m_width = otherPed->m_width;
    m_status = otherPed->m_status;
  }
}
