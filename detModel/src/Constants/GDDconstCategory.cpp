#include "detModel/Constants/GDDconstCategory.h"
#include "detModel/Management/GDDconstantsVisitor.h"

void GDDconstCategory::Accept(GDDconstantsVisitor* v){
  unsigned int i;
  
  v->visitConstCategory(this);
  
  for(i=0; i<consts.size();i++){
    consts[i]->Accept(v);
  }
}

void GDDconstCategory::AcceptNotRec(GDDconstantsVisitor* v){
  v->visitConstCategory(this);
}
