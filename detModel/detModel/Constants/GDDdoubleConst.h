#ifndef GDDDOUBLECONST_H
#define GDDDOUBLECONST_H

#include "detModel/Constants/GDDconst.h"
#include "detModel/Management/GDDconstantsVisitor.h"

class GDDdoubleConst: public GDDconst{
 public:
  GDDdoubleConst():GDDconst(),value(0){}

  virtual void Accept(GDDconstantsVisitor* v){v->visitDoubleConst(this);}
  virtual void AcceptNotRec(GDDconstantsVisitor* v){v->visitDoubleConst(this);}

  double getValue(){return value;}
  void setValue(double pvalue){value=pvalue;}

 private:
 double value;
};

#endif //GDDDOUBLECONST_H
