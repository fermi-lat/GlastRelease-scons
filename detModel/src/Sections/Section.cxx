#include "detModel/Utilities/Purge.h"
#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Sections/Volume.h"
#include "detModel/Sections/Ensemble.h"
#include "detModel/Sections/Shape.h"
#include "detModel/Sections/Section.h"
#include "detModel/Sections/Choice.h"

namespace detModel{

  void Section::Accept(SectionsVisitor* v){
    unsigned int i;
  
    v->visitSection(this);
  
    for(i=0; i<volumes.size();i++){
      volumes[i]->Accept(v);
    }
  }

  void Section::AcceptNotRec(SectionsVisitor* v){
    v->visitSection(this);
  }

  Section::~Section(){
    std::vector<Volume*>::iterator i;
  
    std::vector<Volume*> v = getVolumes();
  
    for(i=v.begin();i!=v.end();i++)
      {
	if(Ensemble* ens = dynamic_cast<Ensemble*>(*i))
	  delete ens;
	else if(Shape* sha = dynamic_cast<Shape*>(*i))
	  delete sha;
	*i = 0;
      }
    purge(choices);
  }

}



