#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Sections/Box.h"
#include "detModel/Sections/BoundingBox.h"

namespace detModel{

  void Box::Accept(SectionsVisitor* v){
    v->visitBox(this);
    if (seg)v->visitSeg(seg);
  }


  void Box::buildBB(){
    getBBox()->setXDim(getX());
    getBBox()->setYDim(getY());
    getBBox()->setZDim(getZ());
  }
  
}
