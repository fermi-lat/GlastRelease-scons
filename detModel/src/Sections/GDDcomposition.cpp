#include "detModel/GDDpurge.h"

#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDcomposition.h"
#include "detModel/Sections/GDDanyPosition.h"
#include "detModel/Sections/GDDposXYZ.h"

GDDcomposition::~GDDcomposition(){
  purge(anyPosition);
}

void GDDcomposition::Accept(GDDsectionsVisitor* v){
    unsigned int i;

    v->visitComposition(this);

    for(i=0; i<anyPosition.size();i++){
      anyPosition[i]->Accept(v);
    }

};


void GDDcomposition::constructBB(){
  unsigned int i;
  double xmin, ymin, zmin, xmax, ymax, zmax;

  xmin = ymin = zmin = 100000;
  xmax = ymax = zmax = -100000;

  for(i=0;i<getPositions().size();i++)
    {
      if (getPositions()[i]->getPosType() == posXYZ)
	{
	  GDDposXYZ* pos = static_cast<GDDposXYZ*>(getPositions()[i]);
      
	  GDDvolume* vol = pos->getVolume();
      
	  if ((pos->getX()+vol->getBBX()/2) > xmax)
	    xmax = (pos->getX()+vol->getBBX()/2);
	  
	  if ((pos->getY()+vol->getBBY()/2) > ymax)
	    ymax = (pos->getY()+vol->getBBY()/2);

	  if ((pos->getZ()+vol->getBBZ()/2) > zmax)
	    zmax = (pos->getZ()+vol->getBBZ()/2);

	  if ((pos->getX()-vol->getBBX()/2) < xmin)
	    xmin = (pos->getX()-vol->getBBX()/2);
	  
	  if ((pos->getY()-vol->getBBY()/2) < ymin)
	    ymin = (pos->getY()-vol->getBBY()/2);

	  if ((pos->getZ()-vol->getBBZ()/2) < zmin)
	    zmin = (pos->getZ()-vol->getBBZ()/2);
	}
      else ///\todo To add other possible positioning
	{
	  xmax = xmin;
	  ymax = ymin;
	  zmax = zmin;
	}
    }//endfor

  bbx = fabs(xmax - xmin);
  bby = fabs(ymax - ymin);
  bbz = fabs(zmax - zmin);
}


