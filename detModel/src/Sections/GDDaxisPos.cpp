#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDaxisPos.h"
#include "detModel/Sections/GDDidField.h"

void GDDaxisPos::Accept(GDDsectionsVisitor* v){
    int i;
    
    v->visitAxisPos(this);
    if(getIdField())
      getIdField()->Accept(v);
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
}
