// $Header$
//   Author: T. Burnett
//  Project: geometry
//

#ifndef GEOMETRY_POINT_H
#define GEOMETRY_POINT_H

#include "geometry/Vector.h"

/** Defines a geometrical point, differing from  Vector only by
  its transformation properties, and the fact that it can be displayed.
  */
class Point: public Hep3Vector , public GeomObject
{
public:
    Point(const Point& p):Hep3Vector(p){}
    Point():Hep3Vector(0,0,0){}
    Point(double x, double y, double z):Hep3Vector(x,y,z){}
    // constructors: default is null
    
    Point& operator=(const Point& v){return (Point&)Hep3Vector::operator=(v);}
    Point& operator+=(const Vector& v){return (Point&)Hep3Vector::operator+=(v);}
    Point& operator-=(const Vector& v){return (Point&)Hep3Vector::operator-=(v);}
    // assignment operators: allow  V = P - P, P = P + V, P = P - V
    
    //(uneeded?) operator Vector()const{return Vector(x(),y(),z());}
    // operator that converts to a vector
    
    GeomObject& transform(const CoordTransform&);
    // transforms like a point
    
    
    const char* nameOf()const{return "Point";}
    void printOn(std::ostream& f =std::cout)const;
    
private:
    // make illegal operations unaccessible that would inherit from Hep3Vector
    Point operator*(double ){return Point();}  // no multiply
    Point operator/(double ){return Point();}  // or divide
    Point operator+(const Point& ){return Point();}; // or add another Point
    
};

// resolve ambiguity bet. GeomObject and Hep3Vector
inline std::ostream& operator<<(std::ostream& out, const Point& v)
{  v.printOn(out); return out;}

// allow  legal addition, subtraction

inline Point operator - (const Point & a, const Vector & b) {
    Point p(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());  return p;}

inline Vector operator - (const Point & a, const Point & b) {
    Vector v(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());  return v;}

inline Point operator + (const Point & a, const Vector & b) {
    Point p(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());  return p;}

inline Point operator + (const Vector & a, const Point & b) {
    Point p(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());  return p;}





#endif

