#ifndef GDDFLOATCONST_H
#define GDDFLOATCONST_H

#include "detModel/Constants/GDDconst.h"
#include "detModel/Management/GDDconstantsVisitor.h"

class GDDfloatConst: public GDDconst{
 public:
  GDDfloatConst():GDDconst(),value(0){}

  virtual void Accept(GDDconstantsVisitor* v){v->visitFloatConst(this);}
  virtual void AcceptNotRec(GDDconstantsVisitor* v){v->visitFloatConst(this);}

  float getValue()const{return value;}
  void setValue(float pvalue){value=pvalue;}

 private:
 float value;
};

#endif //GDDCFLOATCONST_H
