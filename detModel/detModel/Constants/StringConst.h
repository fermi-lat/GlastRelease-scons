#ifndef STRINGCONST_H
#define STRINGCONST_H
#include "detModel/Constants/Const.h"
#include "detModel/Management/ConstantsVisitor.h"

namespace detModel{

  class StringConst: public Const{
  public:
    StringConst():Const(),value(""){}
    
    virtual void Accept(ConstantsVisitor* v){v->visitStringConst(this);}
    virtual void AcceptNotRec(ConstantsVisitor* v){v->visitStringConst(this);}

    std::string getValue()const{return value;}
    void setValue(std::string pvalue){value=pvalue;}
    
  private:
    std::string value;
  };
}

#endif //STRINGCONST_H































