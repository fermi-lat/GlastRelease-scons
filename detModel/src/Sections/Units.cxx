#include "detModel/Sections/Units.h"
#include <string>

namespace detModel{

  void Units::setUnitAngle(std::string s){
    if (s == "deg")
      setUnitAngle(deg);
    else if (s == "rad")
      setUnitAngle(mrad);
  }
  
  void Units::setUnitLength(std::string s){
    if (s == "mm")
      setUnitLength(mm);
    else if (s == "cm")
      setUnitLength(cm);
    else if (s == "m")
      setUnitLength(m);
    
  }

}
