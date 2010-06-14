// $Header$
//

#include "geometry/Hype.h"
#include "geometry/Hyperbolic.h"
#include "geometry/Plane.h"


/*Hype default parameters.
	axis = Vector( 0.0, 0.0, 1.0 );
	half_length = 1.0;
	inner_radius  = 0.9;
	outer_radius = 1.0;
	stereo_angle = 0.1;
*/

Hype::Hype( double hl, double ri, double ro, double sa )
: Volume(4)
{
	if(    ( hl <= 0.0 ) || ( ri <= 0.0 ) || ( ro <= ri )
	     ||( sa <= 0.0 ) || ( sa >=  M_PI/2 )
	  ) FATAL("Illegal Hype parameters");
	Vector axis( 0.0, 0.0, 1.0 );

 	double zt = hl * tan( sa );
	double rio = sqrt( ri * ri  +  zt * zt );
	double roo = sqrt( ro * ro  +  zt * zt );

   //   Create the surfaces: note order!
    addSurface(new Hyperbolic(center(), axis, -ri,-rio,hl));
	addSurface(new Hyperbolic(center(), axis,  ro, roo,hl));
    addSurface(new Plane(center(), Vector(0,0, hl)));
    addSurface(new Plane(center(), Vector(0,0,-hl)));
}

const Vector Hype::axis()const{return ((Hyperbolic&)surface(0)).axis();}

void Hype::printOn( std::ostream& os ) const
{
	Volume::printOn(os);
	os
	   << "axis " << axis()
	   << "\n half_length " << halfLength()
	   << ", inner_radius " << innerRadius()
	   << ", outer_radius " << outerRadius()
	   << ", stereo_angle " << stereoAngle()
	   << "\n"	;
}




