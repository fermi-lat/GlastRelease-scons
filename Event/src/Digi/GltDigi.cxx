// $Header: //

#include "Event/Digi/GltDigi.h"

/** @file 
    @author zachary.fewtrell@nrl.navy.mil

    Implementation file for GltDigi.h
*/

namespace Event {
  unsigned short GltDigi::getCALLETriggerVector() const {
    unsigned short retVal = 0;
    // simple loop through each crystal with high trigger bit & set appropriate bit in return vector
    for (CalTriggerSet::const_iterator it = m_CAL_LO.begin();
         it != m_CAL_LO.end();
         it++) {
      unsigned short tower = it->getTower();
      
      retVal |= 1 << tower;
    }

    return retVal;
  }

  unsigned short GltDigi::getCALHETriggerVector() const {
    unsigned short retVal = 0;
    // simple loop through each crystal with high trigger bit & set appropriate bit in return vector
    for (CalTriggerSet::const_iterator it = m_CAL_HI.begin();
         it != m_CAL_HI.end();
         it++) {
      unsigned short tower = it->getTower();
      
      retVal |= 1 << tower;
    }

    return retVal;
  }

}
