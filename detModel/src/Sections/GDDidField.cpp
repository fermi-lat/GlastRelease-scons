#include <string>
#include "detModel/Sections/GDDidField.h"


void GDDidField::setName(std::string pname){
  name = pname;
}

std::string GDDidField::getName(){
  return name;
}
