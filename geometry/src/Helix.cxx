// $Id$

#include "geometry/Helix.h"

//  static variables for efficiency of position, direction
static double sinSave, cosSave, rhoSave, stepSave, perpSave;


Helix::Helix( const Point& p, const Vector& d,
	      const Vector& a, double r )
	  : Ray( p, d ) ,axis(a.unit()), rho(r)
{  	perp = axis.cross( dir );
	perpSave = perp.magnitude();
	parallel = axis * dir;
	stepSave = sinSave = rhoSave = 0.0;
	cosSave = 1.0;
}

Helix::Helix( const Helix& h ) : Ray( h.pos, h.dir )
{
	axis = h.axis;
	rho = h.rho;
	perp = h.perp;
	parallel = h.parallel;
	arclength = h.arclength;
}

void Helix::printOn( std::ostream& os ) const
{  // printing function using C++ ostream class
	os << "Helix: origin " << Ray::position()
	   << "  dir " << Ray::direction() << "  axis "
	   << axis << "  radius " << rho << "\n";
}
inline void
Helix::updateCache(double step)const
{
   if ( step != stepSave || rho != rhoSave) {
	stepSave = step;
	rhoSave  = rho;
	perpSave = perp.magnitude();
	cosSave = cos( step * perpSave / rho );
	sinSave = sin( step * perpSave / rho );
   }

}

Point Helix::position( double step ) const
{
   updateCache(step);
   double a = parallel * ( step - rho * sinSave / perpSave )
	 ,d = rho * sinSave/perpSave;
   Point p = pos;
   p += d*dir;
   p += a* axis;
   p +=	 ( rho * ( 1.0 - cosSave ) / perpSave ) * perp;
   return p;
 /*  return    pos
		+ d* dir
		+ a * axis
		+ ( rho * ( 1.0 - cosSave ) / perpSave ) * perp
		;
*/
}

Vector
Helix::direction( double step ) const
{
   updateCache(step);
    return    cosSave * dir
		+ parallel * ( 1.0 - cosSave ) * axis
		+ sinSave * perp;
}
GeomObject&
Helix::transform(const CoordTransform& T)
{
   Ray::transform(T);
   axis.transform(T);
   perp.transform(T);
   return *this;
}


