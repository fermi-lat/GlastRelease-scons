#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Sections/AxisMPos.h"
#include "detModel/Sections/BoundingBox.h"

namespace detModel{

  void AxisMPos::Accept(SectionsVisitor* v){
    unsigned int i;
    
    v->visitAxisMPos(this);
    for(i=0;i<getIdFields().size();i++)
      getIdFields()[i]->Accept(v);
  }

  void AxisMPos::translateDisp(double a){
    int i;
  
    for(i=0;i<getNcopy();i++)
      setDisp(i,getDisp(i) + a);
  }

  void AxisMPos::buildBB(){
    int i;

    if (getVolume())
      {
	BoundingBox *b = getVolume()->getBBox();
	if((b->getXDim())*(b->getYDim())*(b->getZDim()) == 0)
	  getVolume()->buildBB();
      
	switch(getAxisDir()){
	case (Stack::xDir):
	  for(i=0;i<getNcopy();i++)
	    addDisp(b->getXDim()/2 + (b->getXDim()+getGap())*i);
	
	  getBBox()->setXDim((b->getXDim()+getGap())*getNcopy() - getGap());
	  getBBox()->setYDim(b->getYDim()); 
	  getBBox()->setZDim(b->getZDim()); 
	  getBBox()->rotate(getRotation(),0,0);
	
	  translateDisp(-(getBBox()->getXDim())/2);
	  break;
	case (Stack::yDir):
	  for(i=0;i<getNcopy();i++)
	    addDisp(b->getYDim()/2 + (b->getYDim()+getGap())*i);

	  getBBox()->setYDim((b->getYDim()+getGap())*getNcopy() - getGap());
	  getBBox()->setXDim(b->getXDim()); 
	  getBBox()->setZDim(b->getZDim()); 
	  getBBox()->rotate(0,getRotation(),0);

	  translateDisp(- getBBox()->getYDim()/2);
	  break;
	case (Stack::zDir):
	  for(i=0;i<getNcopy();i++)
	    addDisp(b->getZDim()/2 + (b->getZDim()+getGap())*i);

	  getBBox()->setZDim((b->getZDim()+getGap())*getNcopy() - getGap());
	  getBBox()->setYDim(b->getYDim()); 
	  getBBox()->setXDim(b->getXDim()); 
	  getBBox()->rotate(0,0,getRotation());

	  translateDisp(- getBBox()->getZDim()/2);
	  break;
	}
      }
  }
}
