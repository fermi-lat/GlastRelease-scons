// $Id$
//


#include "geometry/Cone.h"

#include "geometry/Conic.h"
#include "geometry/Plane.h"
#include <algorithm>
#include <cmath>
inline static double sqr(double x){return x*x;}

static Vector zhat(0,0,1);

Cone::Cone( double length, double ri1, double ro1, double ri2, double ro2 )
: Volume( ((ri1 > 0) && (ri2 >0))? 4:3)
{
  //  Make a Cone with length z, inner radius ri1 and outer radius ro1 at the
  //  -z end, and inner radius ri2, and outer radius ro2 at the +z end.
  //  Note that the following relations must hold for a valid Cone:
  //  z > 0, ri2 > ri1 >= 0, ro2 > ro1 > ri1, and ro2 > ri2.
   if ( ( length <= 0.0 ) || ( ri1 <  0.0 ) || ( ro1 <= ri1 )
                          || ( ri2 <  0.0 ) || ( ro2 <= ri2 )
      )GeometryException("Illegal parameters for Cone");

   double h_z = 0.5*length;
   addSurface(new Conic( center(), zhat, 0.5*(ro1 +ro2), (ro2-ro1)/length) );

   addSurface(new Plane(center(),  h_z*zhat) );
   addSurface(new Plane(center(), -h_z*zhat) );

   if( ri1 > 0 && ri2 >0 ) // add inner surface if not solid
	  addSurface(new Conic( center(), zhat, -0.5*(ri1 +ri2), -(ri2-ri1)/length));

   max_dimension =  sqrt( sqr(length) + sqr(std::max(ro1,ro2)));

}

// This assumes order of surfaces set by code just above. Beware!
double Cone::outerRadius(double z)const
{
    return  ((Conic&)surface(0)).radius(z);
}
double Cone::innerRadius(double z)const
{
    return surfaceCount()>3? -((Conic&)surface(3)).radius(z): 0 ;
}


double Cone::length()const
{
   return 2*((Plane&)surface(1)).offset();
}
Vector Cone::axis()const{return ((Conic&)surface(0)).axis();}

void Cone::printOn( std::ostream& os ) const
{
   double h = 0.5*length();
   Volume::printOn(os);
	os  << "  length " << length()
	    << ", inner radii: " << innerRadius(-h)<< ", " << innerRadius(+h)
	    << ", outer radii: " << outerRadius(-h)<< ", " << outerRadius(+h)
	    << "\n" ;
}


