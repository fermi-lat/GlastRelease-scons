// $Id$
//


#include "geometry/Surface.h"
#include <strstream>

GeomObject&
Surface::transform(const CoordTransform& T )
{
  referenceDirection.transform(T);
  return *this;
}
void
Surface::reverse()
{
    referenceDirection *= -1;
}

// copy constructor uses same reference point, makes a copy of
// direction, then has the subclass do the rest (if it can)
Surface::Surface(const Surface& old)
    :referencePoint(old.referencePoint)
    ,referenceDirection(old.referenceDirection)
{
}

// this to not force all surfaces to copy, but crash if can't
Surface& Surface::copy()const
{
    std::cerr << "The Surface subclass " << nameOf()
	<< " did not implement a copy() method\0";
    return *(Surface*)this; //(just to make compiler happy)
}

