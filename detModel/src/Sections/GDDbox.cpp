#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDbox.h"


void GDDbox::Accept(GDDsectionsVisitor* v){
  v->visitBox(this);
  if (seg)v->visitSeg(seg);
}


void GDDbox::constructBB(){
  setBBX(getX());
  setBBY(getY());
  setBBZ(getZ());
}
