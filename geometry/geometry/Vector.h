// $Header$
//
// Vector.h
//
#ifndef __VECTOR_H
#define __VECTOR_H

#include "CLHEP/Vector/ThreeVector.h"
#include "geometry/GeomObject.h"

/// Subclass of Hep3Vector that has transformations from GeomObject
class Vector : public Hep3Vector , public GeomObject
{
public:
  Vector();
  Vector(double u, double v, double w);
  Vector(double costh, double phi);

  /*explicit*/ Vector(const Hep3Vector& v):Hep3Vector(v){}
  //explicit Vector(const Vector& v){*this = v;}
  // constructors

  int operator==(const Vector& v){
     return x()==v.x() && y()==v.y() && z()==v.z();}

  Vector& operator=(const Vector& v){return (Vector&)Hep3Vector::operator=(v);}
  Vector& operator+=(const Vector& v){return (Vector&)Hep3Vector::operator+=(v);}
  Vector& operator-=(const Vector& v){return (Vector&)Hep3Vector::operator-=(v);}
  Vector operator-()const;

  virtual GeomObject& transform(const CoordTransform& );
 
  GeomObject& rotateX(double a){return GeomObject::rotateX(a);};
  GeomObject& rotateY(double a){return GeomObject::rotateY(a);};
  GeomObject& rotateZ(double a){return GeomObject::rotateZ(a);};
  GeomObject& rotate(double a, const Vector& dir){return GeomObject::rotate(a,dir);}
  // these are to resolve ambiguity with Hep3Vector

  virtual const char* nameOf()const{return "Vector";}


  virtual void printOn( std::ostream& os = std::cout) const;

  double magnitude() const { return mag(); }
  double square() const { return mag2();}

private:
};

inline Vector
Vector::operator-()const{  Vector q(-x(), -y(), -z());  return q;}

// operators that return a vector
inline Vector operator*(double s, const Vector& v){return operator*(s, (const Hep3Vector&) v);}
inline Vector operator*(const Vector& v, double s){return operator*(s, (const Hep3Vector&) v);}

inline Vector operator+(const Vector& a, const Vector& b){return (Vector)((const Hep3Vector&)a+(const Hep3Vector&)b);}
inline Vector operator-(const Vector& a, const Vector& b){return (Vector)((const Hep3Vector&)a-(const Hep3Vector&)b);}

inline std::ostream& operator<<(std::ostream& out, const Vector& v)
{  v.printOn(out); return out;}
// needed to resolve ambiguity bet. GeomObject and Hep3Vector
#endif

