#include "detModel/Utilities/GDDpurge.h"
#include "detModel/Utilities/GDDmatrix.h"
#include "detModel/Sections/GDDstack.h"
#include "detModel/Sections/GDDstackedPos.h"
#include "detModel/Sections/GDDaxisMPos.h"

/// \todo Deal with case different from origin
void GDDstack::buildDisp(){
  double actualDisp = 0;
  
  std::vector<GDDposition*>::const_iterator i;
  GDDboundingBox *b;  
  std::vector<GDDposition*> p = getPositions();
  
  for(i=p.begin();i<p.end();i++)
    {
      b = (*i)->getBBox();
      if(!b->exists())
	(*i)->buildBB();
      
      if(GDDaxisMPos* pos = dynamic_cast<GDDaxisMPos*>(*i))
	{
	  pos->setDispCM(actualDisp + 
			 b->getDirDim((GDDboundingBox::axisDir)getAxisDir())/2 
			 + pos->getGap0());
	  actualDisp = actualDisp + 
	    b->getDirDim((GDDboundingBox::axisDir)getAxisDir()) + pos->getGap0();
	  pos->translateDisp(pos->getDispCM());
	}
    }
}

/// \todo Deal with case different from origin
void GDDstack::buildBB(){
  if(!getBBox()->exists())
    {
      std::vector<GDDposition*>::const_iterator i;
      GDDboundingBox *b;  
      std::vector<GDDposition*> p = getPositions();
      
      buildDisp();
      
      for(i=p.begin();i<p.end();i++)
	{
	  b = (*i)->getBBox();
	  
	  if(GDDaxisMPos* pos = dynamic_cast<GDDaxisMPos*>(*i))
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
	  
	  if(GDDaxisMPos* pos = dynamic_cast<GDDaxisMPos*>(*i))	  
	    {
	      pos->setDispCM(pos->getDispCM() - 
			     getBBox()->getDirDim((GDDboundingBox::axisDir)getAxisDir())/2);
	      pos->translateDisp(-getBBox()->getDirDim((GDDboundingBox::axisDir)getAxisDir())/2);
	    }
	}
    }
}



bool GDDstack::checkPosition(GDDposition* p){
  if (GDDstackedPos* pos = dynamic_cast<GDDstackedPos*>(p))
    {
      pos->setAxisDir(getAxisDir());
      return 1;
    }
  else return 0;
}
