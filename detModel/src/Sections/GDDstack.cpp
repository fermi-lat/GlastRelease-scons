#include "detModel/GDDpurge.h"
#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDstack.h"
#include "detModel/Sections/GDDanyRelativePosition.h"

GDDstack::~GDDstack(){
  purge(anyRelativePosition);
}

void GDDstack::Accept(GDDsectionsVisitor* v){
    unsigned int i;

    v->visitStack(this);

    for(i=0; i<anyRelativePosition.size();i++){
      anyRelativePosition[i]->Accept(v);
    }

};

double GDDstack::getBBX(){
  unsigned int i;
  double res =0;

  switch(getStackType()){
  case sx:
    for(i=0;i<getPositions().size();i++)
      res = res + getPositions()[i]->getBBX();
    break;
  case sy:
  case sz:  
    res = getPositions()[0]->getBBX();
    break;
  }
  return res;
}

double GDDstack::getBBY(){
  unsigned int i;
  double res =0;

  switch(getStackType()){
  case sy:
    for(i=0;i<getPositions().size();i++)
      res = res + getPositions()[i]->getBBY();
    break;
  case sx:
  case sz:  
    res = getPositions()[0]->getBBY();
    break;
  }
  return res;
}

double GDDstack::getBBZ(){
  unsigned int i;
  double res =0;

  switch(getStackType()){
  case sz:
    for(i=0;i<getPositions().size();i++)
      res = res + getPositions()[i]->getBBZ();
    break;
  case sy:
  case sx:  
    res = getPositions()[0]->getBBZ();
    break;
  }
  return res;
}


