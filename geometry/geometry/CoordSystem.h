//     $Id$
// Project: Atlas Geometry
//  Author: Toby Burnett
//  
// This is a coordinate transform that can be transformed: it
// keeps track of the transformations to be able to apply the
// total transformation in either direction.

#ifndef CoordSystem_H
#define CoordSystem_H

#ifdef __GNUG__
#pragma interface
#endif

#include "geometry/GeomObject.h"
#include "geometry/CoordTransform.h"

class CoordSystem : public GeomObject
{

 public:
  CoordSystem();

 
  const CoordTransform& globalToLocal()const;
 // return reference our transform that inverts all the transforms

  const CoordTransform localToGlobal()const;
  // creates (on stack) inverse of the saved globalToLocal xform
 
  GeomObject& transform(const CoordTransform& T);
  //used during construction to establish where we are

 private:
  CoordTransform _globalToLocal;
  // the inverse of the transformation that got our associated
  // Shape to its current position
};

/////////////////////////// Inlines /////////////////////////////////
inline 
CoordSystem::CoordSystem()
: _globalToLocal(CoordTransform::identity)
{}

 
inline GeomObject&
CoordSystem::transform(const CoordTransform& T)
{
   _globalToLocal = T.inverse() * _globalToLocal;
   return *this;
}
inline const CoordTransform&
CoordSystem::globalToLocal()const
{
   return _globalToLocal;
}

inline const CoordTransform 
CoordSystem::localToGlobal()const
{
  return globalToLocal().inverse();
}
#endif      

