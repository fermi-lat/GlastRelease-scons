// $Header$
#ifndef CalibData_HotStrips_h
#define CalibData_HotStrips_h

#include "CalibData/BadStrips.h"

namespace CalibData {
  class HotStrips : virtual public CalibData {
    HotStrips() {
      BadStrips();
      // just invoke base class init 
    }

    // But this doesn't get us anything.  Would a factory pattern
    // somehow  be appropriate? 


    getType() { return HOT; }
  }
}
