// $Header$
#include "CalibData/Anc/AncTaggerGain.h"

namespace CalibData {

  void AncTaggerGain::update(RangeBase* other)  {
    AncTaggerGain* otherGain = dynamic_cast<AncTaggerGain* > (other);
    
    // check that otherGain isn't 0 (dynamic cast worked)
    m_val = otherGain->m_val;
    m_isBad = otherGain->m_isBad;
  }
}
