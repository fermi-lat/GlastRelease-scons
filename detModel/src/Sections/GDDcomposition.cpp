#include "detModel/Utilities/GDDpurge.h"
#include "detModel/Sections/GDDcomposition.h"
#include "detModel/Sections/GDDposition.h"
#include "detModel/Sections/GDDshape.h"
#include "detModel/Sections/GDDposXYZ.h"
#include "detModel/Sections/GDDboundingBox.h"

#include <math.h>

GDDcomposition::~GDDcomposition(){
}

void GDDcomposition::buildBB(){
  std::vector<GDDposition*>::const_iterator i;
  GDDboundingBox *b;  

  std::vector<GDDposition*> p = getPositions();
  
  for(i=p.begin();i<p.end();i++)
    {
      b = (*i)->getBBox();
      if((b->getXDim())*(b->getYDim())*(b->getZDim()) == 0)
	(*i)->buildBB();
      
      if(GDDposXYZ* pos = dynamic_cast<GDDposXYZ*>(*i))
	getBBox()->merge(b,pos->getX(), pos->getY(), pos->getZ());
    }

  b = getEnvelope()->getBBox();
  if((b->getXDim())*(b->getYDim())*(b->getZDim()) == 0)
    getEnvelope()->buildBB();  
  getBBox()->merge(b,0,0,0);
}



