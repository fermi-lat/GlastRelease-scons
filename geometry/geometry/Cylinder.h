// $Id$
//
//
//  Define an cylindrical surface
//
//  Sign of radius parameter determines inside/outside:
//    positive: "inside" is inside the cylinder.

#ifndef __CYLINDER_H
#define __CYLINDER_H

#include "geometry/Surface.h"


class Cylinder: public Surface
{

public:
  Cylinder(const Point& origin, double radius);
  // simple constructor used by Tube: direction on  z-axis
  Cylinder( const Point& origin, const Vector& axis, double radius );
  // constructor: must specify origin point, axis vector, radius
  // where sign of radius indicates inside/or outsice 

  double how_near( const Point& x ) const;
  double distance( const Point& x, const Vector& v, int inout)const;
  // basic specification of this surface

   virtual Vector normal( const Point& x )const;
  // returns normal-vector at Point x, assuming x to be on the surface

  const Vector& axis()const{return Surface::direction();}
  // axis is the surface's basic direction vector 

  const char *nameOf() const { return "Cylinder"; }
  void printOn( std::ostream& os = std::cout ) const;

  double radius()const{return _radius;}
  virtual void setOffset(double offset){_radius=offset;};
  virtual double offset()const{return _radius;};
	

 private:
	double _radius;	// radius of Cylinder (sign indicates inside/outside)
};

#endif

