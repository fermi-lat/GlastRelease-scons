// $Id$
//

// Define a hyperbolic surface of revolution

#include "geometry/Hyperbolic.h"

Hyperbolic::Hyperbolic(const Point& origin , const Vector& axis,
                 double , double , double )
:Surface(origin, axis)
{
   FATAL("not implemented");
}

double
Hyperbolic::how_near( const Point&  ) const
{

  FATAL("not implemented");return 0;
}

double
Hyperbolic::distance( const Point&  ,const Vector& , int) const
{
  FATAL("not implemented");return 0;
}


