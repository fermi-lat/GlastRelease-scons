// $Id$
//
// Define a basic geometrical shape, defined by the inside() and
//  distanceTo... methods
//  This class is abstract

#ifndef __SHAPE_H
#define __SHAPE_H

#include "geometry/Point.h"
#include <float.h> // for FLT_MAX on NT

class Ray;
class Medium;

class Shape : public GeomObject
{
 public:
    virtual double distanceToLeave( const Ray& r, double maxStep  = FLT_MAX )const=0;
    virtual double distanceToEnter( const Ray& r, double maxStep  = FLT_MAX )const=0;
    // return the distance along the Ray to leave (if inside ), or
    // enter (if outside). Can be negative if the origin of the Ray is
    // close to the surface. maxStep is the largest distance required.
    // returns FLT_MAX if no intersection

    virtual int inside ( const Point& x ) const=0;
    // is the point inside the volume?

    void printOn(std::ostream&)const;

    const Point& center()const;
    // access to center data member
    GeomObject& transform(const CoordTransform&);

    virtual Shape& copy()const;
    // non-virtual subclasses figure out how to return a copy

    virtual int getBoundaryIndex(const Point& )const{return -1;}
    // allows client using distanceToLeave to associate an index with a boundary,
    // identified by the point 

    virtual int lastBoundaryIndex()const{return -1;}
    // distanceToLeave sets index of left boundary: client can use

    virtual double getMaxDimension() const { return FLT_MAX; }
    // should be calculated by a subclass which is meant to be the maximum distance 
    // from the center

 protected:
    Shape();
    // default constructor accessible to subclasses

 private:
   Shape(const Shape& ){};
    // copy constructor not allowed

  Point _center;
 // central point
};
inline const Point&
Shape::center()const {return _center;}
#endif



