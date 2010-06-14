//  $Header$
//   Author: T. Burnett
//

#include "geometry/Cons.h"
#include "geometry/Wedge.h"


Cons::Cons( double l, double ri1, double ro1, double ri2, double ro2,
	    double ph1, double ph2 )
: Cone(l,ri1,ro1,ri2,ro2)
{
  addSurface(new Wedge(center(), axis(), ph1, ph2));
}

double Cons::angle1() const
{    int numSurf = surfaceCount();
     Wedge *wedge =   (Wedge *)&surface(numSurf-1);
     const Vector *normal  = &wedge->n1();
     return normal->phi();
}

double Cons::angle2() const
{    int numSurf = surfaceCount();
     Wedge *wedge =   (Wedge *)&surface(numSurf-1);
     const Vector *normal  = &wedge->n2();
     return normal->phi();
}

void Cons::printOn( std::ostream&  ) const
{
}



