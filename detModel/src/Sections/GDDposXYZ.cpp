#include "detModel/Management/GDDsectionsVisitor.h"

#include "detModel/Sections/GDDposXYZ.h"
#include "detModel/Sections/GDDidField.h"

void GDDposXYZ::Accept(GDDsectionsVisitor* v){
    int i;
    
    v->visitPosXYZ(this);
    if(getIdField())
      getIdField()->Accept(v);
    
};


double GDDposXYZ::getBBX(){
  return (getVolume()->getBBX());
}

double GDDposXYZ::getBBY(){
  return (getVolume()->getBBY());
}

double GDDposXYZ::getBBZ(){
  return (getVolume()->getBBZ());
}
