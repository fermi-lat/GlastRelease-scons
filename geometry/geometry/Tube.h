//  $Id$
//
//

#ifndef TUBE_H
#define TUBE_H

#include "geometry/Volume.h"

/// define a Tubular solid
class Tube: public Volume
{

public:
    /// constructor: length, inner, outer radii 
    Tube( double z, double ri, double ro );

    /// resize the Tube with length, inner, outer radii 
    void resize( double z, double ri, double ro );

    const char *nameOf() const { return "Tube"; }
    void printOn( std::ostream& os = std::cout ) const;


    double innerRadius()const;
    double outerRadius()const;
    double length()const;
    // access to dimensions
   
    const Vector& axis()const;	

	
private: 
 
   void calculateSurfaces( double z, double ri, double ro );
};

#endif // TUBE_H

