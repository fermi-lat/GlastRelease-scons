// $Header$
//

//	This volume is a hyperbolic tube designed to be used to describe
//	stereo wires of a drift chamber, for example.
//	The inner and outer surfaces are hyperbolas rotated about the
//	z-axis, while the end surfaces are rings.
//	The user specifies this shape by giving the inner and outer radii
//	(at the center), half length, and stereo angle (in radians).
//	Hype is the same as the GEANT volume of the same name, except that
//	here the stereo angle is given in radians, whereas GEANT gives it
//	in degrees.

#ifndef __HYPE_H
#define __HYPE_H

#include "geometry/Volume.h"


class Hype: public Volume
{
 
public:
    Hype( double z, double ri, double ro, double sa ); 
    // constructor: z is half-length, ri and ro the inner and outer radii,
    // and sa the stereo angle


    const Vector axis()const;
    double halfLength() const { return half_length; }
    double innerRadius() const { return inner_radius; }
    double outerRadius() const { return outer_radius; }
    double stereoAngle() const { return stereo_angle; }

    virtual void calcMaxDimension();


    virtual const char *nameOf() const { return "Hype"; }
    virtual void printOn( std::ostream& os = std::cout ) const;

 private:    // don't really need these, keep them for now          	
    double half_length; 		// half length of Hype
    double inner_radius;		// inner radius of Hype
    double outer_radius;		// outer radius of Hype
    double stereo_angle;		// stereo angle (in radians) of Hype
};

#endif

