// $Id$
//

// Ray subclass that implements a helix
//
#ifndef __HELIX_H
#define __HELIX_H

#include "geometry/Ray.h"

class Helix : public Ray
{
 public:
   Helix( const Point& p, const Vector& d, const Vector& a,
			 double r );
   Helix( const Helix& h );


   Point position( double step ) const;
   Vector direction( double step ) const;
   // position and direction along the helix

   double curvature(double =0) const{ return rho? 1.0/rho : FLT_MAX; }
   // return the local curvature, here just the const. inverse of the radius

   double step() const { return rho? 2.* fabs(rho) : FLT_MAX; }
  // used by Intersection::distance to step along rays 
  // searching for surface crossings

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

