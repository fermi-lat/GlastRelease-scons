#include "detModel/Materials/Element.h"
#include "detModel/Management/MaterialsVisitor.h"

namespace detModel{

  void Element::Accept(MaterialsVisitor* v){
    
    v->visitElement(this);
    
  }

  void Element::AcceptNotRec(MaterialsVisitor* v){
    v->visitElement(this);
  }
  
}
