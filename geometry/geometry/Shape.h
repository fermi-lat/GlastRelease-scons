// $Header$
//

#ifndef GEOMETRY_SHAPE_H
#define GEOMETRY_SHAPE_H

#include "geometry/Point.h"
#include <cfloat> // for FLT_MAX on NT

class Ray;
class Medium;

/**
 Define a basic geometrical shape, defined by the inside() and
  distanceTo... methods
  This class is abstract

*/
class Shape : public GeomObject
{
 public:
     /**
     return the distance along the Ray to leave (if inside ), or
     enter (if outside). Can be negative if the origin of the Ray is
     close to the surface. maxStep is the largest distance required.
     returns FLT_MAX if no intersection
     */
    virtual double distanceToLeave( const Ray& r, double maxStep  = FLT_MAX )const=0;
    virtual double distanceToEnter( const Ray& r, double maxStep  = FLT_MAX )const=0;

    /// is the point inside the volume?
    virtual int inside ( const Point& x ) const=0;

    void printOn(std::ostream&)const;

    /// access to center data member
    const Point& center()const;

    GeomObject& transform(const CoordTransform&);

    /// non-virtual subclasses figure out how to return a copy
    virtual Shape& copy()const;

    /// allows client using distanceToLeave to associate an index with a boundary, identified by the point 
    virtual int getBoundaryIndex(const Point& )const{return -1;}

    /// distanceToLeave sets index of left boundary: client can use
    virtual int lastBoundaryIndex()const{return -1;}

    /// should be calculated by a subclass which is meant to be the maximum distance  from the center
    virtual double getMaxDimension() const { return FLT_MAX; }

    /// how near is it to the nearest surface
    virtual double howNear( const Point& x) const=0;

 protected:
    Shape();
    // default constructor accessible to subclasses

 private:
#ifdef WIN32
   Shape(const Shape& ){};
#else
   Shape(const Shape& ):GeomObject(){}; //gcc wants initializer, msdev chokes
#endif
    // copy constructor not allowed

  Point _center;
 // central point
};
inline const Point&
Shape::center()const {return _center;}
#endif



