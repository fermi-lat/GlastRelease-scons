// $Heading: Cone.h $ 
//   Author: T. Burnett
//
//

#ifndef GEOMETRY_CONE_H
#define GEOMETRY_CONE_H

#include "geometry/Volume.h"

/// Define a right circular conical solid
class Cone: public Volume
{
public:
    /// constructor: specify full length, inner and outer radii at -length/2 inner and outer radii at +length/2
    Cone( double length, double ri1, double ro1, double ri2, double ro2 );

    // Implicit destructor
    // Implicit copy constructor
    // Implicit assignment operator

    Vector axis() const; 
    double length() const; 
    double innerRadius(double) const; 
    double outerRadius(double) const; 
    // access to basic dimensions

    virtual void printOn( std::ostream& os = std::cout ) const;

    virtual const char * nameOf()const{return "Cone";}


};

#endif

