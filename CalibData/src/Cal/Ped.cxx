// $Header$
#include "CalibData/Cal/Ped.h"

namespace CalibData {

  void Ped::update(RangeBase* other)  {
    Ped* otherPed = dynamic_cast<Ped* > (other);
    
    // check that otherPed isn't 0 (dynamic cast worked)
    m_pedAvr = otherPed->m_pedAvr;
    m_pedSig = otherPed->m_pedSig;
  }
}
