// $Id$
//
#include "geometry/Cylinder.h"

inline static double sqr(double x){return x*x;}
#include <float.h>


Cylinder::Cylinder( const Point& origin, const Vector& a, double r )
: Surface( origin, a.unit() )
,_radius(r)
{}
Cylinder::Cylinder(const Point& origin, double r)
:Surface(origin, Vector(0,0,1))
, _radius(r)
{}

double
Cylinder::how_near( const Point& x ) const
{  //  Distance from the point x to the infinite Cylinder.
   //  The distance will be positive if the point is inside the Cylinder,
   //  negative if the point is outside.
	Vector r = x - origin();
	double rho = sqrt( r.mag2() - sqr(r*axis()) );
	return radius()>0?  radius()-rho  : rho+radius() ;
}


Vector
Cylinder::normal( const Point& x )const
{
   Vector r = x-origin();
   Vector rhohat = (r - (r*axis())*axis()).unit(); // unit vector in rho direction
   return radius()>0 ? rhohat : -rhohat;
}
double Cylinder::distance(const Point& x, const Vector& vhat, int inout) const
{
//   solve: | rperp(s) | = R,
//   where rperp(s) = r(s)-ahat*(r(s)*ahat)
//         r(s) = x-origin() + vhat*s

        Vector ahat = axis(),
               r = x - origin(),
	       rperp = r - (ahat*r)*ahat,
	       vperp = vhat-(ahat*vhat)*ahat;

	double a = vperp.mag2(),
	       b = vperp*rperp,
	       c = rperp.mag2() - sqr(radius()),
	       disc = b*b-a*c,
	       s;
        // now want appropriate root of a*s^2 + 2*b*s + c = 0

	if( disc<=0 ) return FLT_MAX;  // misses
	if( radius()<0 ) inout *=-1;	// reverse enter/leave sense
	int leaving = inout==1;
	if( b<0 ) // velocity toward axis
	{  s = (leaving) ? (-b+sqrt(disc))/a :
                          c/(-b+sqrt(disc));
        }
        else // heading away from axis: can only leave
	{  s = (leaving) ? c/(-b-sqrt(disc)) : FLT_MAX;
        }
	return s ;
}

void Cylinder::printOn( std::ostream& os ) const
{	os << "Cylinder surface with origin: " << origin() << "\t"
	   << "radius: " << radius() << "\tand axis " << axis() << "\n";
}

