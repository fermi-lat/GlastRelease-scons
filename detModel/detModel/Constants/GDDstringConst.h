#ifndef GDDSTRINGCONST_H
#define GDDSTRINGCONST_H

#include "detModel/Constants/GDDconst.h"

class GDDstringConst: public GDDconst{
 public:
  GDDstringConst():GDDconst(),value(""){}

  std::string getValue()const{return value;}
  void setValue(std::string pvalue){value=pvalue;}

 private:
 std::string value;
};

#endif //GDDCSTRINGCONST_H































