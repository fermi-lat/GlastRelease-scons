// $Id$
//

#include "geometry/Shape.h"
#include <iostream>

// implementations of Shape member functions

Shape::Shape()
{}


Shape& Shape::copy()const
{

    std::cerr << "The Shape subclass " << nameOf()
	<< " did not implement a copy() method\0";

    return *(Shape*)this; //(just to make compiler happy)

}


void
Shape::printOn(std::ostream& os)const
{
  os << "Shape is a " << nameOf()
	  <<  " centered at: " << center() <<std::endl;
}

GeomObject&
Shape::transform(const CoordTransform& T)
{
  _center.transform(T);
  return *this;
}

