// $Id$
//
//
//      A Trd2 is a list of 2 Rectangles and 4 Trapezoids which bound it.
//      As with GEANT, the Rectangles are in the x-y plane.

#ifndef __TRD2_H
#define __TRD2_H

#include "geometry/Hexahedron.h"

class Trd2: public Hexahedron
{

 public:
    Trd2( double sl, double ll, double sw, double lw, double h );


    double shortLength() const;
    double longLength() const ; 
    double longWidth() const; 
    double height() const; 
    double shortWidth() const; 


    virtual const char *nameOf() const { return "Trd2"; }
    virtual void printOn( std::ostream& os = std::cout ) const;
};

#endif

