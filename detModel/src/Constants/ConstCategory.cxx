#include "detModel/Constants/ConstCategory.h"
#include "detModel/Management/ConstantsVisitor.h"

namespace detModel{
  
  void ConstCategory::Accept(ConstantsVisitor* v){
    unsigned int i;
  
    v->visitConstCategory(this);
    
    for(i=0; i<consts.size();i++){
      consts[i]->Accept(v);
    }
  }

  void ConstCategory::AcceptNotRec(ConstantsVisitor* v){
    v->visitConstCategory(this);
  }
}
