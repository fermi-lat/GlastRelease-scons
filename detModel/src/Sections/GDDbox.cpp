#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDbox.h"
#include "detModel/Sections/GDDboundingBox.h"


void GDDbox::Accept(GDDsectionsVisitor* v){
  v->visitBox(this);
  if (seg)v->visitSeg(seg);
}


void GDDbox::buildBB(){
  getBBox()->setXDim(getX());
  getBBox()->setYDim(getY());
  getBBox()->setZDim(getZ());
}
