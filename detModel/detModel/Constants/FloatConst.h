#ifndef FLOATCONST_H
#define FLOATCONST_H
#include "detModel/Constants/Const.h"
#include "detModel/Management/ConstantsVisitor.h"

namespace detModel{

  class FloatConst: public Const{
  public:
    FloatConst():Const(),value(0){}
    
    virtual void Accept(ConstantsVisitor* v){v->visitFloatConst(this);}
    virtual void AcceptNotRec(ConstantsVisitor* v){v->visitFloatConst(this);}
    
    float getValue()const{return value;}
    void setValue(float pvalue){value=pvalue;}
    
  private:
    float value;
  };
}

#endif //FLOATCONST_H
