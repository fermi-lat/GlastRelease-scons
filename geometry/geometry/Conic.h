//     $Id$
//  Author: T. Burnett
// Project: geometry
//
// Define a conic surface. Inherits origin point and axis direction from Surface

#ifndef GEOMETRY_CONIC_H
#define GEOMETRY_CONIC_H

#include "geometry/Surface.h"

class Conic: public Surface
{
public:

    Conic( const Point& o, const Vector& a, double radius, double slope );
    // specify origin, axis, radius at origin, slope of radius 

    // implicit destructor
    // implicit copy constructor
    // implicit assignment operator

    virtual double how_near( const Point& x ) const;
    virtual double distance( const Point& x ,const Vector& v, int ) const;
    // specify geometry

    virtual Vector normal( const Point& x )const;
    // returns normal-vector at Point x, gradient of howNear function

    const Vector& axis()const{return Surface::direction();}

    double radius(double z=0)const{return m_radius+z*m_slope;}
    // radius at position z from origin along axis

    virtual const char *nameOf() const { return "Conic"; }
    virtual void printOn( std::ostream& os = std::cout ) const;


private:

    double m_slope;    // tangent of opening angle
    double m_radius;    // radius at the point		
    // note that both are negated to reverse sense of surface
};

#endif //GEOMETRY_CONIC_H

