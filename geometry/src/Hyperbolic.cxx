// $Id$
//

// Define a hyperbolic surface of revolution

#include "geometry/Hyperbolic.h"

Hyperbolic::Hyperbolic(const Point& origin , const Vector& axis,
                 double , double , double )
:Surface(origin, axis)
{
   GeometryException("not implemented");
}

double
Hyperbolic::how_near( const Point&  ) const
{

  GeometryException("not implemented");return 0;
}

double
Hyperbolic::distance( const Point&  ,const Vector& , int) const
{
  GeometryException("not implemented");return 0;
}

Vector
Hyperbolic::normal( const Point& ) const
{
  GeometryException("not implemented");
  return Vector();
}

void
Hyperbolic::printOn( std::ostream& ) const
{
  GeometryException("not implemented");
  return;
}


