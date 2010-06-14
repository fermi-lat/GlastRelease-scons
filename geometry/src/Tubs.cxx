// $Header$
//  Author: T. Burnett
//
#include "geometry/Tubs.h"
#include "geometry/Wedge.h"


Tubs::Tubs( double z, double ri, double ro, double phi1, double phi2 )
: Tube(z,ri,ro)
{
  addSurface( new Wedge( center(), axis(),  phi1, phi2 ));
}

const Wedge&
Tubs::wedge()const{
 return (const Wedge&)surface(surfaceCount()-1);}


void Tubs::printOn( std::ostream&  ) const
{
}




