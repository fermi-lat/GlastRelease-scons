#include "detModel/Utilities/Purge.h"
#include "detModel/Utilities/Vector.h"
#include "detModel/Sections/Stack.h"
#include "detModel/Sections/StackedPos.h"
#include "detModel/Sections/AxisMPos.h"

namespace detModel{

  /// \todo Deal with case different from origin
  void Stack::buildDisp(){
    double actualDisp = 0;
  
    std::vector<Position*>::const_iterator i;
    BoundingBox *b;  
    std::vector<Position*> p = getPositions();
  
    for(i=p.begin();i<p.end();i++)
      {
	b = (*i)->getBBox();
	if(!b->exists())
	  (*i)->buildBB();
      
	if(AxisMPos* pos = dynamic_cast<AxisMPos*>(*i))
	  {
	    if(pos->getShift0())
	      {
		pos->setDispCM(actualDisp+pos->getShift0());
		actualDisp = actualDisp+ b->getDirDim((BoundingBox::axisDir)getAxisDir())/2 + 
		  pos->getShift0();			       
	      }
	    else
	      {	    
		pos->setDispCM(actualDisp + 
			       b->getDirDim((BoundingBox::axisDir)getAxisDir())/2 
			       + pos->getGap0());
		actualDisp = actualDisp + 
		  b->getDirDim((BoundingBox::axisDir)getAxisDir()) + pos->getGap0();
	      }
	    pos->translateDisp(pos->getDispCM());
	  }
      }
  }

  /// \todo Deal with case different from origin
  void Stack::buildBB(){
    if(!getBBox()->exists())
      {
	std::vector<Position*>::const_iterator i;
	BoundingBox *b;  
	std::vector<Position*> p = getPositions();
      
	buildDisp();
      
	for(i=p.begin();i<p.end();i++)
	  {
	    b = (*i)->getBBox();
	  
	    if(AxisMPos* pos = dynamic_cast<AxisMPos*>(*i))
	      {
		switch(getAxisDir()){
		case xDir:
		  getBBox()->merge(b,
				   pos->getDispCM()+pos->getDx(), 
				   pos->getDx(), 
				   pos->getDz());
		  getBBox()->setXDim(getBBox()->getXDim()/2);
		  break;	
		case yDir:
		  getBBox()->merge(b,
				   pos->getDx(), 
				   pos->getDispCM()+ pos->getDy(), 
				   pos->getDz());
		  getBBox()->setYDim(getBBox()->getYDim()/2);
		  break;	
		case zDir:
		  getBBox()->merge(b,
				   pos->getDx(), 
				   pos->getDy(),
				   pos->getDispCM()+ pos->getDz());
		  getBBox()->setZDim(getBBox()->getZDim()/2);
		  break;	
		}	
	      }
	  }
	
	for(i=p.begin();i<p.end();i++)
	  {
	    b = (*i)->getBBox();
	  
	    if(AxisMPos* pos = dynamic_cast<AxisMPos*>(*i))	  
	      {
		pos->setDispCM(pos->getDispCM() - 
			       getBBox()->getDirDim((BoundingBox::axisDir)getAxisDir())/2);
		pos->translateDisp(-getBBox()->getDirDim((BoundingBox::axisDir)getAxisDir())/2);
	      }
	  }
      }
  }



  bool Stack::checkPosition(Position* p){
    if (StackedPos* pos = dynamic_cast<StackedPos*>(p))
      {
	pos->setAxisDir(getAxisDir());
	return 1;
      }
    else return 0;
  }
}
