#include "detModel/Utilities/GDDmatrix.h"
#include "detModel/Utilities/GDDpurge.h"
#include "detModel/Sections/GDDboundingBox.h"

#include <iostream>
#include <vector>


void GDDboundingBox::setDirDim(axisDir dir, double a){
  switch(dir){
  case xDir:
    setXDim(a);
    break;
  case yDir:
    setYDim(a);
    break;
  case zDir:
    setZDim(a);
    break;
  }
}
double GDDboundingBox::getDirDim(axisDir dir){
  
  switch(dir){
  case xDir:
    return getXDim();
    break;
  case yDir:
    return getYDim();
    break;
  case zDir:
    return getZDim();
    break;
  }
  return 0;
}


void GDDboundingBox::merge(GDDboundingBox *b, double x, double y, double z){
  double x1min, y1min, z1min, x1max, y1max, z1max;
  double x2min, y2min, z2min, x2max, y2max, z2max;
  double xmin, ymin, zmin, xmax, ymax, zmax;

  x1min = - getXDim()/2;
  y1min = - getYDim()/2;
  z1min = - getZDim()/2;
  
  x1max = getXDim()/2;
  y1max = getYDim()/2;
  z1max = getZDim()/2;

  x2min = (- b->getXDim()/2) + x;
  y2min = (- b->getYDim()/2) + y;
  z2min = (- b->getZDim()/2) + z;
  
  x2max = b->getXDim()/2 + x;
  y2max = b->getYDim()/2 + y;
  z2max = b->getZDim()/2 + z;

  xmin = MIN(x1min,x2min);
  ymin = MIN(y1min,y2min);
  zmin = MIN(z1min,z2min);
  
  xmax = MAX(x1max,x2max);
  ymax = MAX(y1max,y2max);
  zmax = MAX(z1max,z2max);

  setXDim(2*(MAX(fabs(xmax),fabs(xmin))));
  setYDim(2*(MAX(fabs(ymax),fabs(ymin))));
  setZDim(2*(MAX(fabs(zmax),fabs(zmin))));
}

void GDDboundingBox::rotateDir(axisDir dir, double a){
  switch(dir){
  case xDir:
    rotate(a,0,0);
    break;
  case yDir:
    rotate(0,a,0);
    break;
  case zDir:
    rotate(0,0,a);
    break;
  }
}
  
bool GDDboundingBox::exists(){
  if(getXDim()*getYDim()*getZDim() == 0)
    return 0;
  else return 1;
}

void GDDboundingBox::rotate(double xAngle, double yAngle, double zAngle)
{
  unsigned int j;
  std::vector <GDDvector*> points;
  double xmax, ymax, zmax;
 
  xmax = ymax = zmax = -100000;
  
  /// These are 8 points of the bounding box of the positioned volume
  points.push_back(new GDDvector(-getXDim()/2,
				 -getYDim()/2,
				 -getZDim()/2));
  points.push_back(new GDDvector(+getXDim()/2,
				 -getYDim()/2,
				 -getZDim()/2));
  points.push_back(new GDDvector(-getXDim()/2,
				 +getYDim()/2,
				 -getZDim()/2));
  points.push_back(new GDDvector(-getXDim()/2,
				 -getYDim()/2,
				 +getZDim()/2));
  points.push_back(new GDDvector(+getXDim()/2,
				 +getYDim()/2,
				 -getZDim()/2));
  points.push_back(new GDDvector(+getXDim()/2,
				 -getYDim()/2,
				 +getZDim()/2));
  points.push_back(new GDDvector(-getXDim()/2,
				 +getYDim()/2,
				 +getZDim()/2));
  points.push_back(new GDDvector(+getXDim()/2,
				 +getYDim()/2,
				 +getZDim()/2));

  for(j=0;j<8;j++)
      {
	GDDrotate(0,xAngle,points[j]);
	GDDrotate(1,yAngle,points[j]);	  
	GDDrotate(2,zAngle,points[j]);
      }
  
  for(j=0;j<8;j++)
    {
      if (points[j]->x > xmax)
	xmax = points[j]->x;
      
      if (points[j]->y > ymax)
	ymax = points[j]->y;
      
      if (points[j]->z > zmax)
	zmax = points[j]->z;
    }

  setXDim(2*xmax);
  setYDim(2*ymax);
  setZDim(2*zmax);
  
  purge(points);
}
