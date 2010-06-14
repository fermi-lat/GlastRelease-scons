//  $Id$
//
//

#ifndef __SPHERE_H
#define __SPHERE_H

#include "geometry/Surface.h"

/// Define a spherical surface. Inherit origin point from Surface
class Sphere: public Surface
{
    
public:
    // specify origin, radius 
    Sphere( const Point& o,  double radius );
    
    
    virtual double how_near( const Point& x ) const;
    virtual double distance( const Point& x ,const Vector& v, int ) const;

    virtual Vector normal( const Point& x )const;
    // returns normal-vector at Point x, assuming x to be "on" the surface
    // specify geometry
    
    
    double radius()const{return _radius;}
    // radius
    
    virtual const char *nameOf() const { return "Sphere"; }
    virtual void printOn( std::ostream& os = std::cout ) const;
    
    
private:
    
    double _radius;
    // radius 		
};

#endif

