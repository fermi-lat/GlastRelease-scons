#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDaxisPos.h"
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Sections/GDDboundingBox.h"


void GDDaxisPos::Accept(GDDsectionsVisitor* v){
  unsigned int i;
    
  v->visitAxisPos(this);
  for(i=0;i<getIdFields().size();i++)
    getIdFields()[i]->Accept(v);
}


void GDDaxisPos::buildBB(){
  if (getVolume())
    {
      GDDboundingBox *b = getVolume()->getBBox();
      if((b->getXDim())*(b->getYDim())*(b->getZDim()) == 0)
	getVolume()->buildBB();
      
      getBBox()->setXDim(b->getXDim()); 
      getBBox()->setYDim(b->getYDim()); 
      getBBox()->setZDim(b->getZDim()); 
      
      switch(getAxisDir()){
      case xDir:
	getBBox()->rotate(getRotation(),0,0);
	break;
      case yDir:
	getBBox()->rotate(0,getRotation(),0);
	break;
      case zDir:
	getBBox()->rotate(0,0,getRotation());
	break;
      }
    }
}

