#ifndef GDDSTRINGCONST_H
#define GDDSTRINGCONST_H

#include "detModel/Constants/GDDconst.h"
#include "detModel/Management/GDDconstantsVisitor.h"

class GDDstringConst: public GDDconst{
 public:
  GDDstringConst():GDDconst(),value(""){}

  virtual void Accept(GDDconstantsVisitor* v){v->visitStringConst(this);}
  virtual void AcceptNotRec(GDDconstantsVisitor* v){v->visitStringConst(this);}

  std::string getValue()const{return value;}
  void setValue(std::string pvalue){value=pvalue;}

 private:
 std::string value;
};

#endif //GDDCSTRINGCONST_H































