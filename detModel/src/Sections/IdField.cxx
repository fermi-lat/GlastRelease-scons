#include <string>
#include "detModel/Sections/IdField.h"

namespace detModel{

  void IdField::setName(std::string pname){
    name = pname;
  }

  std::string IdField::getName()const{
    return name;
  }
}
