#ifndef GDDSTRINGCONST_H
#define GDDSTRINGCONST_H

#include "detModel/Constants/GDDconst.h"

class GDDstringConst: public GDDconst{
 public:
  GDDstringConst():GDDconst(),value(""){}

  string getValue()const{return value;}
  void setValue(string pvalue){value=pvalue;}

 private:
 string value;
};

#endif //GDDCSTRINGCONST_H































