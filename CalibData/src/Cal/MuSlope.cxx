// $Header$
#include "CalibData/Cal/MuSlope.h"

namespace CalibData {

  void MuSlope::update(RangeBase* other)  {
    MuSlope* otherMuSlope = dynamic_cast<MuSlope* > (other);
    
    // check that otherMuSlope isn't 0 (dynamic cast worked)
    m_muSlope = otherMuSlope->m_muSlope;
  }
}
