// CoordTransform.cxx,v 2.1 1995/05/17 23:52:58 burnett Exp
//


#include "geometry/CoordTransform.h"
#include "geometry/Point.h"

// allocate a standard identity transform
const CoordTransform CoordTransform::identity;

CoordTransform
CoordTransform::inverse()const
{
  return CoordTransform(R.inverse(),-(R.inverse()*X));

}
CoordTransform
CoordTransform::operator*(const CoordTransform& T) const
{
  return CoordTransform((T.R  * R) , /*(Vector&)*/ (T.R*X) + T.X);
}

void
CoordTransform::transformCoord(Point& p)const {

  Hep3Vector q(p);
  q *= R;
  //  p = (Point&)(R*p);
  p = Point( q.x()+X.x(), q.y()+X.y(), q.z()+X.z() );

}

