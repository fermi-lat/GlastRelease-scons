#include "detModel/Utilities/GDDpurge.h"
#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Sections/GDDensamble.h"
#include "detModel/Sections/GDDshape.h"
#include "detModel/Sections/GDDchoice.h"
#include "detModel/Sections/GDDsection.h"

void GDDsection::Accept(GDDsectionsVisitor* v){
  unsigned int i;
  
  
  v->visitSection(this);
  
  for(i=0; i<volumes.size();i++){
    volumes[i]->Accept(v);
  }
}

void GDDsection::AcceptNotRec(GDDsectionsVisitor* v){
  v->visitSection(this);
}

GDDsection::~GDDsection(){
  std::vector<GDDvolume*>::iterator i;
  
  std::vector<GDDvolume*> v = getVolumes();
  
  for(i=v.begin();i!=v.end();i++)
    {
      if(GDDensamble* ens = dynamic_cast<GDDensamble*>(*i))
	delete ens;
      else if(GDDshape* sha = dynamic_cast<GDDshape*>(*i))
	delete sha;
      *i = 0;
    }
  purge(choices);
}





