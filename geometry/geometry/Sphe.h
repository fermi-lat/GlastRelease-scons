// $Header$
//  Author: T. Burnett
//

//	C++ header file for the  Sphe class which is a segment of a
//	sphere defined by inner and outer radii and azimuthal and polar
//	angle limits.
//	It is not clear if this Sphe is the same as the GEANT volume of the
//	same name.

#ifndef __SPHE_H
#define __SPHE_H

#include "geometry/Volume.h"

class Sphe: public Volume
{
    
public:
    Sphe( double ri,  double ro,
        double th1=0, double th2=M_PI,
        double ph1=0, double ph2=2*M_PI );
    // constructor: inner, outer radii;
    // limits in theta, phi (defaults: full solid angle)
    
    
    double innerRadius() const{return inner_radius;}
    double outerRadius() const{return outer_radius;}
    double theta1() const{return theta_1;}
    double theta2() const{return theta_2;}
    double phi1() const{return phi_1;}
    double phi2() const{return phi_2;}
    // access to construction parameters
    
    int fullAzimuth() const;
    int fullPolar() const;
    // flags for convenience
    
    virtual const char *nameOf() const { return "Sphe"; }
    virtual void printOn( std::ostream& os = std::cout ) const;
private: // keep these for now
    double inner_radius;	// inner radius of Sphe, can be zero
    double outer_radius;	// outer radius of Sphe, must be greater
				// than inner radius
    double theta_1;		// lower polar angle limit of Sphe
				// in radians, allowed range is
				// 0 <= theta_1 < PI
    double theta_2;		// upper polar angle limit of Sphe
				// in radians, allowed range is
				// theta_1 < theta_2 <= theta_1 + PI
    double phi_1;   	// lower azimuthal angle limit of Sphe
				// in radians, allowed range is
				// 0 <= phi_1 < 2*PI
    double phi_2;		// upper azimuthal angle limit of Sphe
				// in radians, allowed range is
				// phi_1 < phi_2 <= phi_1 + 2*PI
    
    
    int full_azimuth;	// flag set to 1 if the azimuthal angular
				// limits cover the full 2 PI 	
    int full_polar;		// flag set to 1 if the polar angular
				// limits cover the full PI
				// These flags are set by the constructor
				// and resize functions. 	
    
};

#endif

