// $Header$
#ifndef GEOMETRY_RAY_H
#define GEOMETRY_RAY_H

#include "geometry/Point.h"
#include <cfloat> // needed on NT for FLT_MAX

/// GeomObject that is a base for directed rays, implementing a straight line
class Ray : public GeomObject
{
public:
  Ray( const Point& p, const Vector& d );
  Ray( const Ray& r );
  // constructor

  /// return position along ray
  virtual Point position( double s ) const ;
  /// position at start of ray
  const Point& position() const {return pos;}

  /// return direction along ray
  virtual Vector direction( double  ) const{ return dir; }
  /// direction at start of ray
  const Vector& direction() const {return dir;}

  /// transform the Ray to a new system
  GeomObject& transform( const CoordTransform & );

  /// allows polymorphic test to check local curvature
  virtual double curvature(double =0) const { return 0.0; }

  /// used by Intersection::distance to step along rays searching for surface crossings
  virtual double step() const {return FLT_MAX;}

  double getArcLength() const { return arclength; }
  void setArcLength( double s ) { arclength = s; }
  void setFlag( int k ) {flag = k;}
  int  getFlag() const {return flag;} 
  void addArcLength( double ds ) { arclength += ds; }
  // access to internal position

  /// useful for vertexing between this Ray and another
  virtual double distanceOfClosestApproach(const Ray& r) const;
  virtual Vector vectorOfClosestApproach(const Ray& r) const;

  const char *nameOf() const { return "Ray"; }
  virtual void printOn( std::ostream& os = std::cout ) const;


 protected:
  Point pos;			// starting point
  Vector dir;   		// direction unit vector
  float arclength;	// length along the Ray to be drawn
  int flag;                     // a flag to be used by clients

};

#endif

