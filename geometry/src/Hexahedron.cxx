// $Id$
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "geometry/Hexahedron.h"

inline static double sqr(double x){return x*x;}

Hexahedron::Hexahedron()
: Volume(6)
{}

Point
Hexahedron::vertex(int n)const
{ // map to the plane to use:
  int ix=(n&1), iy= ((n/2)&1)+2, iz=((n/4)&1)+4;
  Point dummy;
  Vector a(plane(ix).normal(dummy)),
         b(plane(iy).normal(dummy)),
         c(plane(iz).normal(dummy));
  double l=plane(ix).offset(), w=plane(iy).offset(), h=plane(iz).offset();

  double ax=a.x(), ay=a.y(), az=a.z(),
        bx=b.x(), by=b.y(), bz=b.z(),
	cx=c.x(), cy=c.y(), cz=c.z();
  // invert the 3x3 matrix by hand
  double Mxx = by*cz-bz*cy,
    Myx =-bx*cz+bz*cx,
	Mzx = bx*cy-by*cx,

	Mxy =-ay*cz+az*cy,
	Myy = ax*cz-az*cx,
	Mzy =-ax*cy+ay*cx,

	Mxz = ay*bz-az*by,
	Myz =-ax*bz+az*bx,
	Mzz = ax*by-ay*bx,

	D = ax*by*cz-ax*bz*cy-bx*ay*cz+bx*az*cy+cx*ay*bz-cx*az*by;

   if( abs(D)<1e-10 )	  // probably bad code if this happens
      return Point(0,0,0);
   Vector v( (Mxx*l + Mxy*w + Mxz*h)/D,
             (Myx*l + Myy*w + Myz*h)/D,
             (Mzx*l + Mzy*w + Mzz*h)/D);
// Next line commented out & replaced because gcc tries to make a Point from a Vector.
//   return center()+v;
     Point p = center();
     return p+v;
}

