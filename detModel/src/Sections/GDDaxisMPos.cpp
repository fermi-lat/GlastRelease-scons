#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Sections/GDDaxisMPos.h"
#include "detModel/Sections/GDDidField.h"

void GDDaxisMPos::Accept(GDDsectionsVisitor* v){
  unsigned int i;
    
  v->visitAxisMPos(this);
  for(i=0;i<getIdFields().size();i++)
    getIdFields()[i]->Accept(v);
}


double GDDaxisMPos::getBBX(){
  switch(getPosDir()){
  case(sx):
     return (getVolume()->getBBX() + getGap())*getNcopy() -getGap(); 
    break;
  case(sy):
  case(sz):
     return (getVolume()->getBBX());
    break;
  }
  return 0;
}

double GDDaxisMPos::getBBY(){
  switch(getPosDir()){
  case(sy):
    return (getVolume()->getBBY()+getGap())*getNcopy() - getGap(); 
    break;
  case(sx):
  case(sz):
     return (getVolume()->getBBY());
    break;
  }
  return 0;
}

double GDDaxisMPos::getBBZ(){
  switch(getPosDir()){
  case(sz):
    return (getVolume()->getBBZ() + getGap())*getNcopy() - getGap(); 
    break;
  case(sx):
  case(sy):
    return (getVolume()->getBBZ());
    break;
  }
  return 0;
}


