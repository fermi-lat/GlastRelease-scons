// $Header$
//

#include "geometry/Point.h"
#include "geometry/CoordTransform.h"



GeomObject&
Point::transform(const CoordTransform& T)
{
	T.transformCoord(*this);
    return *this;
}




void
Point::printOn(std::ostream& cout)const
{
	cout << "Point";
	cout << '(' << x() << ',' << y() << ',' << z() << ')';
}


