#include "detModel/Sections/GDDensemble.h"
#include "detModel/Sections/GDDposition.h"
#include "detModel/Management/GDDsectionsVisitor.h"

#include <vector>

void GDDensemble::Accept(GDDsectionsVisitor* v){
  std::vector<GDDposition*>::const_iterator i;
  std::vector<GDDposition*> p = getPositions();
  
  v->visitEnsemble(this);
  for(i=p.begin(); i<p.end();i++){
    (*i)->Accept(v);
  }
};
