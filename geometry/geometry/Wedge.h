// $Header$
//

#ifndef __WEDGE_H
#define __WEDGE_H

#include "geometry/Surface.h"
#include "geometry/Plane.h"

/// Define wedge-like surface formed by two intesecting planes
class Wedge : public Surface
{
 public:
    Wedge(const Point& c,const Vector& axis, double phi1, double phi2);
    // constructor: inside is from phi1 to phi2

    virtual double how_near( const Point& x )const;
    // signed function giving minimum distance to the surface: surface
    // is actually defined by zero.

    virtual double distance( const Point& x ,const Vector& v, int inout) const;
    virtual Vector normal( const Point& x )const;
    // returns normal-vector at Point x, assuming x to be "on" the surface

    const Vector& axis()const{return Surface::direction();}
    const Vector& n1()const; //return the normal to plane 1
    const Vector& n2()const; // return the normal to plane 2
 
    const char * nameOf()const{return "Wedge";}

 private:
    Vector _n1;
    Vector _n2;
    double dphi;
};
#endif

