// $Id$
//

#include "geometry/Vector.h"

#include "geometry/CoordTransform.h"

Vector::Vector(double costh, double phi)
: Hep3Vector(sqrt(1-costh*costh)*cos(phi)
	    ,sqrt(1-costh*costh)*sin(phi)
	    ,costh)
	    {}
// convenient constructor of unit vector in given direction?

Vector::Vector()
  : Hep3Vector(){}

Vector::Vector(double u, double v, double w)
  : Hep3Vector(u,v,w){}

GeomObject&
Vector::transform(const CoordTransform& T)
{
  T.transformVector(*this);
  return *this;
}


void
Vector::printOn( std::ostream& os ) const
{
	os << '(' << x() << ',' << y() << ',' << z() << ')';
}



