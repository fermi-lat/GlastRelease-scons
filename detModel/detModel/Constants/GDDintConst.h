#ifndef GDDINTCONST_H
#define GDDINTCONST_H

#include "detModel/Constants/GDDconst.h"

class GDDintConst: public GDDconst{
 public:
  GDDintConst():GDDconst(),value(0){;}

  int getValue()const{return value;}
  void setValue(int pvalue){value=pvalue;}

 private:
 int value;
};

#endif //GDDCINTCONST_H
