// $Header$
#include "CalibData/Acd/AcdPed.h"

namespace CalibData {

  void AcdPed::update(RangeBase* other)  {
    AcdPed* otherPed = dynamic_cast<AcdPed* > (other);
    
    // check that otherPed isn't 0 (dynamic cast worked)
    m_ped = otherPed->m_ped;
    m_sig = otherPed->m_sig;
  }
}
