// $Id$
//
#include "geometry/Trd1.h"

static Vector xhat(1,0,0),yhat(0,1,0), zhat(0,0,1);

Trd1::Trd1( double sl, double ll, double w, double h )
: Hexahedron()
{
  //  Make a Trd1 with short_length (sl) and long_length (ll) along x-axis,
  //  width (w) along y-axis, and height (h) along z-axis.
  //  As with GEANT, the 2 Trapezoids are in the x-z plane.

  //  Require the long_length to be longer than or equal to the short_length
  if ( ( sl <= 0.0 ) || ( ll <= 0.0 ) || ( w <= 0.0 ) || ( h <= 0.0 )
         || ( ll < sl )
     ) FATAL("Illegal parameters in Trd1");

  max_dimension = 0.5 * sqrt( ll*ll + w*w + h*h );

  // +x, -x planes: note z-component
  double a = atan(0.5*(ll-sl)/h ),
         s = sin(a),
         c = cos(a),
         d = c* 0.25*(ll+sl) ;// perp distance to planes

  addPlane( d * Vector( c, 0, -s));
  addPlane( d * Vector(-c, 0, -s));

  // +y, -y planes:
  addPlane(  0.5*w * yhat);
  addPlane( -0.5*w * yhat);

  // +z, -z planes:
  addPlane(  0.5*h * zhat);
  addPlane( -0.5*h * zhat);
}


double Trd1::shortLength() const{return 0;}	 // fix when someone complains
double Trd1::longLength()  const{return 0;}
double Trd1::width()  const{return 0;}
double Trd1::height() const{return 0;}

void Trd1::printOn( std::ostream& os ) const
{
    Volume::printOn(os);
   	 os    << "short_length "  << shortLength()
           << ", long_length " << longLength()
           << "\n width "      << width()
           << ", height "      << height()
           << "\n";
}




