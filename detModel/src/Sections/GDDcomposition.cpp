#include "detModel/Utilities/GDDpurge.h"
#include "detModel/Utilities/GDDmatrix.h"
#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDcomposition.h"
#include "detModel/Sections/GDDanyPosition.h"
#include "detModel/Sections/GDDposXYZ.h"

#include <math.h>

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
  unsigned int i, j;
  double xmin, ymin, zmin, xmax, ymax, zmax;

  
  xmin = ymin = zmin = 100000;
  xmax = ymax = zmax = -100000;

  for(i=0;i<getPositions().size();i++)
    {
      if (getPositions()[i]->getPosType() == posXYZ)
	{
	  std::vector <GDDvector*> points;
	  GDDposXYZ* pos = static_cast<GDDposXYZ*>(getPositions()[i]);
      	  GDDvolume* vol = pos->getVolume();
	  
	  /// These are 8 points of the bounding box of the positioned volume
	  points.push_back(new GDDvector(-vol->getBBX()/2,
					-vol->getBBY()/2,
					-vol->getBBZ()/2));

	  points.push_back(new GDDvector(+vol->getBBX()/2,
					-vol->getBBY()/2,
					-vol->getBBZ()/2));

	  points.push_back(new GDDvector(-vol->getBBX()/2,
					+vol->getBBY()/2,
					-vol->getBBZ()/2));

	  points.push_back(new GDDvector(-vol->getBBX()/2,
					-vol->getBBY()/2,
					+vol->getBBZ()/2));

	  points.push_back(new GDDvector(+vol->getBBX()/2,
					+vol->getBBY()/2,
					-vol->getBBZ()/2));
	  points.push_back(new GDDvector(+vol->getBBX()/2,
					-vol->getBBY()/2,
					+vol->getBBZ()/2));
	  points.push_back(new GDDvector(-vol->getBBX()/2,
					+vol->getBBY()/2,
					+vol->getBBZ()/2));

	  points.push_back(new GDDvector(+vol->getBBX()/2,
					+vol->getBBY()/2,
					+vol->getBBZ()/2));

	  GDDvector* t = new GDDvector(pos->getX(),pos->getY(),pos->getZ());
	  
	  
	  if(pos->getZrot() != 0)
	    for(j=0;j<8;j++)
	      GDDrotate(2,pos->getZrot(),points[j]);
	  if(pos->getYrot() != 0)
	    for(j=0;j<8;j++)
	      GDDrotate(1,pos->getYrot(),points[j]);	  
	  if(pos->getXrot() != 0)
	    for(j=0;j<8;j++)
	      GDDrotate(0,pos->getXrot(),points[j]);


	  for(j=0;j<8;j++)
	    GDDtranslate(t,points[j]);


	  for(j=0;j<8;j++)
	    {
	      if (points[j]->x > xmax)
		xmax = points[j]->x;

	      if (points[j]->y > ymax)
		ymax = points[j]->y;

	      if (points[j]->z > zmax)
		zmax = points[j]->z;

	      if (points[j]->x < xmin)
		xmin = points[j]->x;

	      if (points[j]->y < ymin)
		ymin = points[j]->y;

	      if (points[j]->z < zmin)
		zmin = points[j]->z;
	    }
	  purge(points);
	}///\todo To add other possible positioning
    }//endfor

  setOX(xmin);
  setOY(ymin);
  setOZ(zmin);
  setBBX(fabs(xmax - xmin));
  setBBY(fabs(ymax - ymin));
  setBBZ(fabs(zmax - zmin));
}


