#ifndef GDDFLOATCONST_H
#define GDDFLOATCONST_H

#include "detModel/Constants/GDDconst.h"

class GDDfloatConst: public GDDconst{
 public:
  GDDfloatConst():GDDconst(),value(0){}

  float getValue()const{return value;}
  void setValue(float pvalue){value=pvalue;}

 private:
 float value;
};

#endif //GDDCFLOATCONST_H
