// $Header$
//
#ifndef GEOMETRY_HELIX_H
#define GEOMETRY_HELIX_H

#include "geometry/Ray.h"

/// Ray subclass that implements a helix
class Helix : public Ray
{
 public:
   Helix( const Point& p, const Vector& d, const Vector& a,
			 double r );
   Helix( const Helix& h );


   /// position along the Helix
   Point position( double step ) const;
   /// direction at the given distance along the Helix
   Vector direction( double step ) const;

   /// return the local curvature, here just the const. inverse of the radius
   double curvature(double =0) const{ return rho? 1.0/rho : FLT_MAX; }


  /// used by Intersection::distance to step along rays searching for surface crossings
   double step() const { return rho? 2.* fabs(rho) : FLT_MAX; }

   GeomObject& transform( const CoordTransform &);

   const char *nameOf() const { return "Helix"; }
   void printOn( std::ostream& os = std::cout ) const;

 private:
   Vector axis;		// helix axis direction (unit vector)
   double rho; 		// helix radius, sign significant
   Vector perp;  	// perpendicular direction
   double parallel;	// component along axis (for efficiency)
   void updateCache(double)const;
};

#endif

