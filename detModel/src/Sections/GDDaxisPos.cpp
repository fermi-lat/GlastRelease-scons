#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDaxisPos.h"
#include "detModel/Sections/GDDidField.h"

void GDDaxisPos::Accept(GDDsectionsVisitor* v){
  unsigned int i;
    
  v->visitAxisPos(this);
  for(i=0;i<getIdFields().size();i++)
    getIdFields()[i]->Accept(v);
}

double GDDaxisPos::getBBX(){
  return getVolume()->getBBX(); 
}

double GDDaxisPos::getBBY(){
  return getVolume()->getBBY(); 
}

double GDDaxisPos::getBBZ(){
  return getVolume()->getBBZ(); 
}

