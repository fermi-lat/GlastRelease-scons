#include "detModel/Sections/Ensemble.h"
#include "detModel/Sections/Position.h"
#include "detModel/Management/SectionsVisitor.h"

#include <vector>

namespace detModel{

  void Ensemble::Accept(SectionsVisitor* v){
    std::vector<Position*>::const_iterator i;
    std::vector<Position*> p = getPositions();
    
    v->visitEnsemble(this);
    for(i=p.begin(); i<p.end();i++){
      (*i)->Accept(v);
    }
  };

}
