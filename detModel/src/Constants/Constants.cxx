#include "detModel/Constants/Constants.h"
#include "detModel/Management/ConstantsVisitor.h"

namespace detModel{

  void Constants::Accept(ConstantsVisitor* v){
    unsigned int i;
    
    v->visitConstants(this);
    
    for(i=0; i<constCategories.size();i++){
      constCategories[i]->Accept(v);
    }
  }

  void Constants::AcceptNotRec(ConstantsVisitor* v){
    v->visitConstants(this);
  }
  
}
