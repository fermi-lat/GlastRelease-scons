// $Id$

#include "geometry/Ray.h"

Ray::Ray( const Point& p, const Vector& d )
: pos( p ), dir(d.unit()) , arclength(0), flag(0)
{
}

Ray::Ray( const Ray& r )
:pos(r.pos),GeomObject(), dir(r.dir), arclength(r.arclength)
{  // copy constructor
}


Point
Ray::position(double s) const
{
//   return  pos + (s * dir);
  Point p = pos;
  return p + s*dir;
}
//should define operators to do this


void
Ray::printOn( std::ostream& os ) const
{  // printing function using C++ ostream class
	os << "Ray: origin = " << pos << ", direction = "
	   << dir << "\n";
}

GeomObject&
Ray::transform( const CoordTransform & T)
{
  pos.transform(T);
  dir.transform(T);
  return *this;
}

double
Ray::distanceOfClosestApproach( const Ray& r ) const
{
   Vector perp = direction(0.).cross(r.direction(0.));
   double d = (perp*(r.position(0.)-pos)) / perp.mag();
   return d;
}

Vector
Ray::vectorOfClosestApproach( const Ray& r ) const
{
   Vector perp = direction(0.).cross(r.direction(0.));
   double d = (perp*(r.position(0.)-pos)) / perp.mag();
   return d*perp;
}



