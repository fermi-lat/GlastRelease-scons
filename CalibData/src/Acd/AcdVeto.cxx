// $Header$
#include "CalibData/Acd/AcdVeto.h"

namespace CalibData {

  void AcdVeto::update(RangeBase* other)  {
    AcdVeto* otherVeto = dynamic_cast<AcdVeto* > (other);
    
    // check that otherVeto isn't 0 (dynamic cast worked)
    m_veto = otherVeto->m_veto;
    m_width = otherVeto->m_width;
    m_status = otherVeto->m_status;
  }
}
