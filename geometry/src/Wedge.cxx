// Wedge.cxx,v 2.1 1995/05/17 23:54:13 burnett Exp
//

#include "geometry/Wedge.h"
#include <cfloat>

Wedge::Wedge(const Point& center, const Vector& axis,  double phi1, double phi2)
:Surface(center, axis )
, _n1( sin(phi1),-cos(phi1), 0)
, _n2(-sin(phi2), cos(phi2), 0)
, dphi(phi2-phi1)
{
    // constructor: inside is from phi1 to phi2  ?

    // phi1 is angle from x-axis to first plane
    // phi2 is angle from x-axis to second plane
	if ( ( phi2 < phi1 ) || ( phi2 >= ( phi1 + 2*M_PI ) ) )
	  GeometryException("Bad parameters for a Wedge");

}

const Vector&
Wedge::n1() const {   return _n1; }

const Vector&
Wedge::n2() const {   return _n2; }

double
Wedge::how_near( const Point& x ) const
{

    // Initializes planes to test
    Plane plane1(origin(), n1(), 0),
          plane2(origin(), n2(), 0);

    double d1 = plane1.how_near(x),
           d2 = plane2.how_near(x),
           line; // distance from the intersection line

    Vector r = x-origin(),
           rperp = r - (r*axis())*axis();

    line = rperp.mag();

    int type = (dphi <= M_PI),
        pl1  = (d1 > 0),
        pl2  = (d2 > 0);

   if (type)
   {
       // four possibilities: Point is in front of both, behind both,
       // or in front of one and not the other(2x)
       if (pl1 && pl2)
       {
           return (d1>d2)? d2:d1;
       }
       else if (!pl1 && !pl2)
       {
           return -line;
       }
       else if (!pl1 && pl2)
       {
          // Now it gets tricky...

          // distance will not always be the shortest distance of line and plane,
          // because only a certain region of the plane exists.
          // To decide which to use, we must determine if the point x
          // is over the existant part of the plane, or the non-existant part.
          // The point will be over the existant part if some vector
          // if it lies in front of a plane A that is perpendicular to
          // plane1, intersects at the intersection line, and faces the existant
          // part of plane1.

          // step 1: find an arbitrary vector that is perpendicular to
          // both the normal of 1 and the axis

          Vector dir = n1().cross(axis());

          // Now we need to see if this vector is in the right direction.
          // This vector should be the sAme direction as the other normal
          // If not, reverse it
          // NOte: If dphi > M_PI, then it is reverse of the other normal
          double n = dir*n2();
          if (n<0)
             dir = -dir;

          // Now dir is prepared.

          n = r*dir; // Is r the same direction as n?

          return (n>0)? d1:-line; // If it is, then the distance t the plane is correct.
                                // because it lies over the existant part of the plane
       }
       else if (pl1 && !pl2)
       {
          // Same as the procedure above
          Vector dir = n2().cross(axis());

          double n = dir*n1();

          if (n<0)
             dir = -dir;

          n = r*dir;

          return (n>0)? d2:-line;
       }
   }
   else // this means dphi > M_PI
   {
       if (pl1 && pl2)
       {
          return line;
       }
       else if (!pl1 && !pl2)
       {
          return (d1>d2)? d2:d1;
       }
       else if (!pl1 && pl2)
       {
          Vector dir = n2().cross(axis());

          double n = dir*n1();

          if (n>0)
             dir = -dir;

          n = r*dir;

          return (n>0)? d2:line; // Remember: return value must be positive because it is out
       }
       else if (pl1 && !pl2)
       {
          Vector dir = n1().cross(axis());

          double n = dir*n2();

          if (n>0)
             dir = -dir;


          n = r*dir;

          return (n>0)? d1:line;
       }

    }
   return(FLT_MAX);
}


double
Wedge::distance(const Point& x, const Vector& v, int inout)const
{

  Plane p[]={Plane(origin(), n1(), 0),
             Plane(origin(), n2(), 0)};

  double d[] ={ p[0].distance(x,v,inout),
                p[1].distance(x,v,inout) };

  if( d[0]==FLT_MAX && d[1]==FLT_MAX )
     return FLT_MAX;

  unsigned solution = d[0]<d[1]?0:1;
  double dist = d[solution];

  if( dphi< M_PI )
  {
    Point hit(x);
    hit += dist*v;
    if( p[1-solution].inside(hit) )
      return dist;
    else
      return FLT_MAX;
  }
  return FLT_MAX;
}

Vector
Wedge::normal( const Point& x )const
{

    Plane p[]={Plane(origin(), n1(), 0),
             Plane(origin(), n2(), 0)};

    double d[] ={ p[0].how_near(x),
                p[1].how_near(x) };

    return fabs(d[0])<fabs(d[1]) ? n1() : n2() ;
}







