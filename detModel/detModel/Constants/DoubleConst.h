#ifndef DOUBLECONST_H
#define DOUBLECONST_H
#include "detModel/Constants/Const.h"
#include "detModel/Management/ConstantsVisitor.h"

namespace detModel{

  class DoubleConst: public Const{
  public:
    DoubleConst():Const(),value(0){}
    
    virtual void Accept(ConstantsVisitor* v){v->visitDoubleConst(this);}
    virtual void AcceptNotRec(ConstantsVisitor* v){v->visitDoubleConst(this);}

    double getValue(){return value;}
    void setValue(double pvalue){value=pvalue;}
    
  private:
    double value;
  };
}

#endif //DOUBLECONST_H
