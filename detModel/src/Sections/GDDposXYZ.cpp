#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Sections/GDDposXYZ.h"
#include "detModel/Sections/GDDboundingBox.h"

void GDDposXYZ::Accept(GDDsectionsVisitor* v){
  unsigned int i;
  v->visitPosXYZ(this);
  for(i=0;i<getIdFields().size();i++)
    getIdFields()[i]->Accept(v);
}

void GDDposXYZ::buildBB(){
  if (getVolume())
    {
      GDDboundingBox *b = getVolume()->getBBox();
      if((b->getXDim())*(b->getYDim())*(b->getZDim()) == 0)
	getVolume()->buildBB();
      
      getBBox()->setXDim(b->getXDim()); 
      getBBox()->setYDim(b->getYDim()); 
      getBBox()->setZDim(b->getZDim()); 
      
      getBBox()->rotate(getXRot(),getYRot(),getZRot());
    }
}

