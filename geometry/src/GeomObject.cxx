// $Id$
//
//

#include "geometry/GeomObject.h"
#include "geometry/CoordTransform.h"

#include <typeinfo>

const char* GeomObject::nameOf()const
{
	const type_info& t1 = typeid(*this);
	return t1.name();
}

GeomObject&
GeomObject::move(const Vector& v)
{
  return transform(CoordTransform(v));
}
GeomObject& GeomObject::moveX(double s){ return move(Vector(s,0,0));}
GeomObject& GeomObject::moveY(double s){ return move(Vector(0,s,0));}
GeomObject& GeomObject::moveZ(double s){ return move(Vector(0,0,s));}

GeomObject&
GeomObject::rotateX(double theta)
{
    HepRotation R;
    return transform(R.rotateX(theta));
}


GeomObject &
GeomObject::rotateY(double theta)
{
    HepRotation R;
    return transform(R.rotateY(theta));
}

GeomObject &
GeomObject::rotateZ(double theta)
{
    HepRotation R;
    return transform(R.rotateZ(theta));
}

GeomObject &
GeomObject::rotate(double angle, const Vector & axis)
{
    HepRotation R;
    return transform(R.rotate(angle, axis));
}
void GeomObject::printOn(std::ostream& )const{};





