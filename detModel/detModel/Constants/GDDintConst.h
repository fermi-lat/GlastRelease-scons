#ifndef GDDINTCONST_H
#define GDDINTCONST_H

#include "detModel/Constants/GDDconst.h"
#include "detModel/Management/GDDconstantsVisitor.h"

class GDDintConst: public GDDconst{
 public:
  GDDintConst():GDDconst(),value(0){;}

  virtual void Accept(GDDconstantsVisitor* v){v->visitIntConst(this);}
  virtual void AcceptNotRec(GDDconstantsVisitor* v){v->visitIntConst(this);}

  int getValue()const{return value;}
  void setValue(int pvalue){value=pvalue;}

 private:
 int value;
};

#endif //GDDCINTCONST_H
