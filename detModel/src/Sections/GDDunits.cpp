#include "detModel/Sections/GDDunits.h"
#include <string>

void GDDunits::setUnitAngle(std::string s){
  if (s == "deg")
    setUnitAngle(deg);
  else if (s == "rad")
    setUnitAngle(mrad);
}

void GDDunits::setUnitLength(std::string s){
  if (s == "mm")
    setUnitLength(mm);
  else if (s == "cm")
    setUnitLength(cm);
  else if (s == "m")
    setUnitLength(m);

}
