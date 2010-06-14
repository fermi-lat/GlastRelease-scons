// $Id$
//

#include "geometry/Trd2.h"


Trd2::Trd2( double sl, double ll, double sw, double lw, double h )
: Hexahedron()
{
//  Make a Trd2 with short_length (sl) and long_length (ll) along x-axis,
//  short_width (sw) and long_width along y-axis,
//  and height (h) along z-axis.
    if ( ( sl <= 0.0 ) || ( ll <= 0.0 ) ||
             ( sw <= 0.0 ) || ( lw <= 0.0 ) || ( h <= 0.0 )
			|| ( ll < sl )
			|| ( lw < sw )
	   ) FATAL("Illegal parameters for Trd2");
	max_dimension = 0.5 * sqrt( ll*ll  + lw*lw + h*h );

    double ldif = 0.5 *(ll - sl),
           wdif = 0.5 *(lw - sw),
           lsum = 0.25*(lw + sw),
           wsum = 0.25*(lw + sw),
           lh = sqrt( h*h + ldif * ldif ),
           wh = sqrt( h*h + wdif * wdif );

    // Planes at +x, -x
	double c = h/lh,
	       s = ldif/lh,
	       d = lsum*c; // distance to x-planes
	addPlane( d*Vector( c,0,-s) );
	addPlane( d*Vector(-c,0,-s) );

	// Planes at +y, -y
    c = h/wh; s = wdif/wh; d = wsum*c;
	addPlane( d* Vector(0, c,-s) );
	addPlane( d* Vector(0,-c,-s) );

    // Planes at +z, -z
	addPlane( Vector(0,0,  0.5*h) );
	addPlane( Vector(0,0, -0.5*h) );

}

double Trd2::shortLength() const{return 0;}	 // fix when someone complains
double Trd2::longLength()  const{return 0;}
double Trd2::shortWidth()  const{return 0;}
double Trd2::longWidth()  const{return 0;}
double Trd2::height() const{return 0;}


void Trd2::printOn( std::ostream& os ) const
{
   Volume::printOn(os);
    os  << " short_length " << shortLength()
        << ",long_length "  << longLength()
        << "\n  short_width " << shortWidth()
        << ", long_width " << longWidth()
        << ", and height " << height()
        << "\n";
}




