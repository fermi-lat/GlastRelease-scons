// $Id$
//
//
//  Trd1 is the same as the GEANT volume of the same name, but note that
//  GEANT uses half-lengths as parameters, whereas we full lengths.
//      A Trd1 has a List of 4 Rectangles and 2 Trapezoids which bound it.
//      As with GEANT, the Trapezoids are in the x-z plane.
//

#ifndef __TRD1_H
#define __TRD1_H

#include "geometry/Hexahedron.h"


class Trd1: public Hexahedron
{


public:
    Trd1( double shortLength, double longLength, double width, double height );
	// constructor:
	//  shortLength  and longLength  along x-axis,
    //  width  along y-axis, and height along z-axis.


    double shortLength() const;
    double longLength() const ; 
    double width() const; 
    double height() const; 

    virtual const char *nameOf() const { return "Trd1"; }
    virtual void printOn( std::ostream& os = std::cout ) const;


};

#endif

