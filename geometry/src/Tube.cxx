// $Id$
//


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "geometry/Tube.h"

#include "geometry/Cylinder.h"
#include "geometry/Plane.h"

inline static double sqr(double x){return x*x;}

static Vector zhat(0,0,1);

Tube::Tube( double length, double innerRadius, double outerRadius )
: Volume( (innerRadius>0)? 4 : 3)
{
    calculateSurfaces(length, innerRadius, outerRadius);
}

void Tube::resize( double length, double innerRadius, double outerRadius )
{
    deleteSurfaces();
    calculateSurfaces(length, innerRadius, outerRadius);
}


void Tube::calculateSurfaces( double length, double innerRadius, double outerRadius )
{
   if(    innerRadius<0 || innerRadius>outerRadius
       || outerRadius<0 || length <=0
      )GeometryException("Illegal parameters for Tube::calculateSurfaces");

    // don't reorder these without looking below!
	addSurface(new Cylinder(center(),  outerRadius));
    addSurface(new Plane(center(), 0.5*length*zhat));
    addSurface(new Plane(center(),-0.5*length*zhat));
    if( innerRadius !=0)
	addSurface(new Cylinder(center(), -innerRadius));  // must be last!	
	
    max_dimension =  sqrt( sqr(length) + sqr(outerRadius));
}
// These assume order of surfaces set by code just above. Beware!

double Tube::outerRadius()const{ return ((Cylinder&)surface(0)).radius();}

double Tube::innerRadius()const{
    return surfaceCount()>3? -((Cylinder&)surface(3)).radius(): 0 ;}

double Tube::length()const
{  return  2* ((Plane&)surface(1)).offset();
}
const Vector& Tube::axis()const{return ((Cylinder&)surface(0)).axis();}



void Tube::printOn( std::ostream& os ) const
{  //  printing function, uses C++ ostream class object
	Volume::printOn(os);
	os  << "  length " << length()
	    << ", inner, outer radii: "
	    << innerRadius()<< ", " << outerRadius() << "\n" ;
}





