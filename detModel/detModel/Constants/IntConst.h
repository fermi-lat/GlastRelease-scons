#ifndef INTCONST_H
#define INTCONST_H
#include "detModel/Constants/Const.h"
#include "detModel/Management/ConstantsVisitor.h"

namespace detModel{

  class IntConst: public Const{
  public:
    IntConst():Const(),value(0){;}
    
    virtual void Accept(ConstantsVisitor* v){v->visitIntConst(this);}
    virtual void AcceptNotRec(ConstantsVisitor* v){v->visitIntConst(this);}
    
    int getValue()const{return value;}
    void setValue(int pvalue){value=pvalue;}
    
  private:
    int value;
  };
}

#endif //INTCONST_H
