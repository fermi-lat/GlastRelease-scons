#include "detModel/Constants/GDDconstants.h"
#include "detModel/Management/GDDconstantsVisitor.h"

void GDDconstants::Accept(GDDconstantsVisitor* v){
  unsigned int i;
  
  v->visitConstants(this);
  
  for(i=0; i<constCategories.size();i++){
    constCategories[i]->Accept(v);
  }
}

void GDDconstants::AcceptNotRec(GDDconstantsVisitor* v){
  v->visitConstants(this);
}
