// $Header$
//  Author: T. Burnett
//
#include "geometry/Sphe.h"

#include "geometry/Sphere.h"


Sphe::Sphe( double ri, double ro,
	double th1, double th2,
	double ph1, double ph2 )
	: Volume(3)  //THB: temporary to allow cap to function
{
//  Make a Sphe with inner radius ri, outer radius ro, lower polar limit
//  th1, upper polar limit th2, lower azimuthal limit ph1, and upper azimuthal
//  limit ph2.

	if ( ( ri < 0.0 ) || ( ro <= ri ) || ( th1 < 0.0 ) || ( th1 >= M_PI )
	  || ( th2 <= th1 ) || ( th2 > ( th1 + M_PI ) ) || ( ph1 < 0.0 )
	  || ( ph1 >= 2*M_PI ) || ( ph2 <= ph1 ) || ( ph2 > ( ph1 + 2*M_PI ) ) )
	{
		GeometryException("Dimensions of Sphe are invalid");
	}
	inner_radius  = ri;
	outer_radius  = ro;
	theta_1 = th1;
	theta_2 = th2;
	phi_1 = ph1;
	phi_2 = ph2;


//  Set the flags to indicate how many surfaces will be needed to describe
//  this Sphe.
	full_polar =  fabs( theta_2 - theta_1 - M_PI ) <= FLT_EPSILON ;
	full_azimuth =  fabs( phi_2 - phi_1 - 2*M_PI ) <= FLT_EPSILON ;

	if( !full_polar || !full_azimuth )
	   GeometryException("Partial sphe not yet implemented, sorry");

	addSurface(new Sphere(center(), ro) );
	if( ri>0 )
	   addSurface(new Sphere(center(), -ri) );

}


void Sphe::printOn( std::ostream& os ) const
{
  Shape::printOn(os);
  os   << "\t radius range:      " << innerRadius() << " - " << outerRadius() << "\n"
       << "\t polar angle range: " << theta1() << " - " << theta2() << " radians\n"
       << "\t azimuthal range:   " << phi1() << " - " << phi2() << " radians\n";

}

