// $Id$
//
// CoordTransform.h
//
//  Define basic coordinate transformations 
//
//   Vector: v' = R * v
//   Point:  p' = R * p + X
//
// where R is a 3x3 rotation matrix and X a Vector.

#ifndef __COORDTRANSFORM_H
#define __COORDTRANSFORM_H

#include "CLHEP/Vector/Rotation.h"
#include "geometry/Vector.h"
class Point;

typedef HepRotation Rotation;

class CoordTransform
{
  public:
   CoordTransform(){};
   CoordTransform(const Vector& t):X(t){};
   CoordTransform(const Rotation& _R):R(_R){};
   CoordTransform(const Rotation& _R, const Vector& t):X(t),R(_R){}
   // constructors: specify either or both rotation matrix, translation

   CoordTransform inverse()const;
   // return inverse transformation

   CoordTransform operator*(const CoordTransform& T) const;
   // product of two transformations

   void transformVector(Vector& v)const{ v*=R;}
   void transformCoord(Point& p)const;

   static const CoordTransform identity;
   // special default transform 
  private:

   Vector X;   // translation  vector
   Rotation R; // rotation  matrix
};

// simple little inline to get around conversion problems with some compilers
inline Vector operator*(const Rotation& R, const Vector& p)
{
  Hep3Vector q(p);
  return Vector(R*q);
}
#endif

