// $Id$
// GeomObject that is a base for directed rays, implementing a straight line
#ifndef __RAY_H
#define __RAY_H

#include "geometry/Point.h"
#include <float.h> // needed on NT for FLT_MAX

class Ray : public GeomObject
{
public:
  Ray( const Point& p, const Vector& d );
  Ray( const Ray& r );
  // constructor

  virtual Point position( double s ) const ;
  const Point& position() const {return pos;}
  // return position at start, or along ray

  virtual Vector direction( double  ) const{ return dir; }
  const Vector& direction() const {return dir;}
  // return direction at start, or along ray

  GeomObject& transform( const CoordTransform & );
  // transform the Ray to a new system

  virtual double curvature(double =0) const { return 0.0; }
  // allows polymorphic test to check local curvature

  virtual double step() const {return FLT_MAX;}
  // used by Intersection::distance to step along rays 
  // searching for surface crossings

  double getArcLength() const { return arclength; }
  void setArcLength( double s ) { arclength = s; }
  void addArcLength( double ds ) { arclength += ds; }
  // access to internal position

  virtual double distanceOfClosestApproach(const Ray& r) const;
  virtual Vector vectorOfClosestApproach(const Ray& r) const;
  // useful for vertexing between this Ray and another

  const char *nameOf() const { return "Ray"; }
  virtual void printOn( std::ostream& os = std::cout ) const;


 protected:
  Point pos;			// starting point
  Vector dir;   		// direction unit vector
  float arclength;	// length along the Ray to be drawn

};

#endif

