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
  switch(getPosDir()){
  case(sx):
    return getVolume()->getBBX() + getGap(); 
    break;
  case(sy):
  case(sz):
    return (getVolume()->getBBX());
    break;
  }
  return 0;
}

double GDDaxisPos::getBBY(){
  switch(getPosDir()){
  case(sy):
    return getVolume()->getBBY() + getGap(); 
    break;
  case(sx):
  case(sz):
    return (getVolume()->getBBY());
    break;
  }
  return 0;
}

double GDDaxisPos::getBBZ(){
  switch(getPosDir()){
  case(sz):
    return getVolume()->getBBZ() + getGap(); 
    break;
  case(sx):
  case(sy):
    return (getVolume()->getBBZ());
    break;
  }
  return 0;
}
