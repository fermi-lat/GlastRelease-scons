#include "detModel/Sections/GDDensamble.h"
#include "detModel/Sections/GDDposition.h"
#include "detModel/Management/GDDsectionsVisitor.h"

#include <vector>

void GDDensamble::Accept(GDDsectionsVisitor* v){
  std::vector<GDDposition*>::const_iterator i;
  std::vector<GDDposition*> p = getPositions();
  
  v->visitEnsamble(this);
  for(i=p.begin(); i<p.end();i++){
    (*i)->Accept(v);
  }
};
