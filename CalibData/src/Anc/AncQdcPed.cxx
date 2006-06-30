// $Header$
#include "CalibData/Anc/AncQdcPed.h"

namespace CalibData {

  void AncQdcPed::update(RangeBase* other)  {
    AncQdcPed* otherPed = dynamic_cast<AncQdcPed* > (other);
    
    // check that otherPed isn't 0 (dynamic cast worked)
    m_val = otherPed->m_val;
    m_rms = otherPed->m_rms;
    m_isBad = otherPed->m_isBad;
    m_device = otherPed->m_device;
  }
}
