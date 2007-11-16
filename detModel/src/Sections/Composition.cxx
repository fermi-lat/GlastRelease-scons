#include "detModel/Utilities/Purge.h"
#include "detModel/Sections/Composition.h"
#include "detModel/Sections/Shape.h"
#include "detModel/Sections/SinglePos.h"
#include "detModel/Sections/PosXYZ.h"
#include "detModel/Sections/BoundingBox.h"

namespace detModel{

  Composition::~Composition(){
  }

  void Composition::buildBB(){
    std::vector<Position*>::const_iterator i;
    BoundingBox *b;  

    std::vector<Position*> p = getPositions();
  
    for(i=p.begin();i<p.end();i++)
      {
	b = (*i)->getBBox();
	if((b->getXDim())*(b->getYDim())*(b->getZDim()) == 0)
	  (*i)->buildBB();
      
	if(PosXYZ* pos = dynamic_cast<PosXYZ*>(*i))
	  getBBox()->merge(b,pos->getX(), pos->getY(), pos->getZ());
      }

    b = getEnvelope()->getBBox();
    if((b->getXDim())*(b->getYDim())*(b->getZDim()) == 0)
      getEnvelope()->buildBB();  
    getBBox()->merge(b,0,0,0);
  }

  bool Composition::checkPosition(Position* p){
    if (dynamic_cast<SinglePos*>(p))
      return 1;
    else return 0;
  }

}
