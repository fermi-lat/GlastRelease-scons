#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Sections/Volume.h"
#include "detModel/Sections/PosXYZ.h"
#include "detModel/Sections/BoundingBox.h"

namespace detModel{

  void PosXYZ::Accept(SectionsVisitor* v){
    unsigned int i;
    v->visitPosXYZ(this);
    for(i=0;i<getIdFields().size();i++)
      getIdFields()[i]->Accept(v);
  }
  
void PosXYZ::buildBB(){
  if (getVolume())
    {
      BoundingBox *b = getVolume()->getBBox();
      if((b->getXDim())*(b->getYDim())*(b->getZDim()) == 0)
	getVolume()->buildBB();
      
      getBBox()->setXDim(b->getXDim()); 
      getBBox()->setYDim(b->getYDim()); 
      getBBox()->setZDim(b->getZDim()); 
      
      getBBox()->rotate(getXRot(),getYRot(),getZRot());
    }
}
}

