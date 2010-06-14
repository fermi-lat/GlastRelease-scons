// $Id$
//
//
#include "geometry/Sphere.h"

#include "CLHEP/config/TemplateFunctions.h"
#include <cfloat>

Sphere::Sphere( const Point& origin, double radius)
: Surface(origin, Vector(0,0,0) )
, _radius(radius)
{}


void Sphere::printOn( std::ostream&  ) const
{
}

double
Sphere::how_near( const Point& x ) const
{  //  Distance from the point x to the nearest point on the Sphere.
   //  The distance will be positive if the point is inside the Sphere,
   //  negative if the point is outside.
 	Vector r = x - origin();
	double rabs = r.mag();

	return radius()>0?  radius()-rabs  : rabs+radius() ;
 }


Vector
Sphere::normal( const Point& x )const
{
    // returns normal-vector at Point x, assuming x to be "on" the surface,
    // and not coinciding with the origin!

    if(radius()>0 ) return (x - origin()).unit() ;
    else            return (origin() - x).unit() ;
}
double
Sphere::distance(const Point& x, const Vector& v, int inout)const
{
   // note very similar logic to Cylinder::distance
   Vector r = x-origin();
   double b = r*v,
          c= r.mag2()-sqr(radius()),
	  disc = b*b-c,
	  s;

   // now want appropriate root of s^2 + 2*b*s + c = 0

   if( disc<=0 ) return FLT_MAX;  // misses
   if( radius()<0 ) inout *=-1;	  // reverse enter/leave sense
   int leaving = inout==1;
   if( b<0 ) // velocity toward axis
   {  s = (leaving) ? -b+sqrt(disc) :
                         c/(-b+sqrt(disc));
   }
   else // heading away from axis: can only leave
   {  s = (leaving) ? c/(-b-sqrt(disc)) : FLT_MAX;
   }
  return s>0? s : FLT_MAX;

}

