#include "detModel/Utilities/GDDpurge.h"
#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDvolume.h"
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
  purge(volumes);
  purge(choices);
}






