// $Header$
//  Author: T. Burnett
//
// Define a Tube with a cut

#ifndef __TUBS_H
#define __TUBS_H

#include "geometry/Tube.h"
#include "geometry/Wedge.h"

class Tubs: public Tube
{

public:
    Tubs( double z, double ri, double ro, double ph1, double ph2 );

    virtual const char *nameOf() const { return "Tubs"; }
    virtual void printOn( std::ostream& os = std::cout ) const;


    const Wedge& wedge() const;
    // access to the Wedge that defines the Tubs
};

#endif

