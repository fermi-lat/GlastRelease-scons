#ifndef GDDDOUBLECONST_H
#define GDDDOUBLECONST_H

#include "detModel/Constants/GDDconst.h"

class GDDdoubleConst: public GDDconst{
 public:
  GDDdoubleConst():GDDconst(),value(0){}

  double getValue()const{return value;}
  void setValue(double pvalue){value=pvalue;}

 private:
 double value;
};

#endif //GDDDOUBLECONST_H
