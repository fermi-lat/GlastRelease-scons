#include "detModel/GDDpurge.h"

#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDcomposition.h"
#include "detModel/Sections/GDDanyPosition.h"

GDDcomposition::~GDDcomposition(){
  purge(anyPosition);
}

void GDDcomposition::Accept(GDDsectionsVisitor* v){
    int i;

    v->visitComposition(this);

    for(i=0; i<anyPosition.size();i++){
      anyPosition[i]->Accept(v);
    }

};


double GDDcomposition::getBBX(){
  return getEnvelope()->getBBX();
}

double GDDcomposition::getBBY(){
  return getEnvelope()->getBBY();  
}

double GDDcomposition::getBBZ(){
  return getEnvelope()->getBBZ();
}


