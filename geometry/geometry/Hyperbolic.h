// $Id$
//

#ifndef HYPERBOLIC_H
#define HYPERBOLIC_H

#include "geometry/Surface.h"

class Hyperbolic : public Surface
{
  public:

    Hyperbolic(const Point& origin, const Vector& axis, double ri, double ro, double h_z);

    virtual double how_near( const Point& x ) const;
    virtual double distance( const Point& x ,const Vector& v,int) const;
    virtual Vector normal( const Point& x )const;
    // returns normal-vector at Point x, assuming x to be "on" the surface

    // specify geometry

    const Vector& axis()const{return Surface::direction();}

    double radius(double z=0)const;
    // radius at position z from origin

    virtual const char *nameOf() const { return "Hyperbolic"; }
    virtual void printOn( std::ostream& os = std::cout ) const;




};
#endif

