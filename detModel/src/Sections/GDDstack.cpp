#include "detModel/Utilities/GDDpurge.h"
#include "detModel/Utilities/GDDmatrix.h"
#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDstack.h"
#include "detModel/Sections/GDDanyRelativePosition.h"

#include <math.h>

GDDstack::~GDDstack(){
  purge(anyRelativePosition);
}

void GDDstack::Accept(GDDsectionsVisitor* v){
    unsigned int i;

    v->visitStack(this);

    for(i=0; i<anyRelativePosition.size();i++){
      anyRelativePosition[i]->Accept(v);
           }

};

void GDDstack::constructBB(){
  unsigned int i, j;
  double xmin, ymin, zmin, xmax, ymax, zmax;
  std::vector <GDDvector*> points;
  double trans,step = 0;

  xmin = ymin = zmin = 100000;
  xmax = ymax = zmax = -100000;


  points.push_back(new GDDvector(0.0,0.0,0.0));
  for(i=0;i<getPositions().size();i++)
    {
      step=0;  
      /// These are 8 points of the bounding box of the positioned volume
      points.push_back(new GDDvector(-getPositions()[i]->getBBX()/2,
				     -getPositions()[i]->getBBY()/2,
				     -getPositions()[i]->getBBZ()/2));
      
      points.push_back(new GDDvector(+getPositions()[i]->getBBX()/2,
				     -getPositions()[i]->getBBY()/2,
				     -getPositions()[i]->getBBZ()/2));
      
      points.push_back(new GDDvector(-getPositions()[i]->getBBX()/2,
				     +getPositions()[i]->getBBY()/2,
				     -getPositions()[i]->getBBZ()/2));
      
      points.push_back(new GDDvector(-getPositions()[i]->getBBX()/2,
				     -getPositions()[i]->getBBY()/2,
				     +getPositions()[i]->getBBZ()/2));
      
      points.push_back(new GDDvector(+getPositions()[i]->getBBX()/2,
				     +getPositions()[i]->getBBY()/2,
				     -getPositions()[i]->getBBZ()/2));
      points.push_back(new GDDvector(+getPositions()[i]->getBBX()/2,
				     -getPositions()[i]->getBBY()/2,
				     +getPositions()[i]->getBBZ()/2));
      points.push_back(new GDDvector(-getPositions()[i]->getBBX()/2,
				     +getPositions()[i]->getBBY()/2,
				     +getPositions()[i]->getBBZ()/2));
      
      points.push_back(new GDDvector(+getPositions()[i]->getBBX()/2,
				     +getPositions()[i]->getBBY()/2,
				     +getPositions()[i]->getBBZ()/2));
      
      switch(getStackType()){
      case sx:
	if (getPositions()[i]->getShift()==0){
	  step = getPositions()[i]->getGap() + getPositions()[i]->getBBX()/2 + getPositions()[i]->getDx();
	  if (i>0)step+=getPositions()[i-1]->getBBX()/2;
	}
	else
	  step=getPositions()[i]->getShift()+getPositions()[i]->getDx();
	trans+=step;
	for(j=points.size()-8;j<points.size();j++){
	  if(getPositions()[i]->getRotation() != 0)
	    GDDrotate(0,getPositions()[i]->getRotation(),points[j]);
	  GDDtranslate(new GDDvector(trans,getPositions()[i]->getDy(),getPositions()[i]->getDz()), points[j]);
	}	
	break;
      case sy:
	if (getPositions()[i]->getShift()==0){
	  step = getPositions()[i]->getGap() + getPositions()[i]->getBBY()/2 + getPositions()[i]->getDy();
	  if (i>0)step+=getPositions()[i-1]->getBBY()/2;
	}
	else
	  step=getPositions()[i]->getShift()+getPositions()[i]->getDy();
	trans+=step;
	for(j=points.size()-8;j<points.size();j++){
	  if(getPositions()[i]->getRotation() != 0)
	    GDDrotate(1,getPositions()[i]->getRotation(),points[j]);
	  GDDtranslate(new GDDvector(getPositions()[i]->getDx(),trans,getPositions()[i]->getDz()),points[j]);
	}	
	break;
      case sz:
	if (getPositions()[i]->getShift()==0){
	  step = getPositions()[i]->getGap() + getPositions()[i]->getBBZ()/2 + getPositions()[i]->getDz();
	  if (i>0)step+=getPositions()[i-1]->getBBZ()/2;
	}
	else
	  step=getPositions()[i]->getShift()+getPositions()[i]->getDz();
	trans+=step;
	for(j=points.size()-8;j<points.size();j++){
	  if(getPositions()[i]->getRotation() != 0)
	    GDDrotate(2,getPositions()[i]->getRotation(),points[j]);
	  GDDtranslate(new GDDvector(getPositions()[i]->getDx(),getPositions()[i]->getDy(),trans),points[j]);
	}	
	break;
      }//end switch
    }//endfor
  

  for(j=0;j<points.size();j++)
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
  setOX(xmin);
  setOY(ymin);
  setOZ(zmin);
  setBBX(fabs(xmax - xmin));
  setBBY(fabs(ymax - ymin));
  setBBZ(fabs(zmax - zmin));
}

