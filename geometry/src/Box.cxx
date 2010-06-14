// $Id$
//
//
// Implement a box: this is it, folks!


#include "geometry/Box.h"

static Vector _xhat(1,0,0), _yhat(0,1,0), _zhat(0,0,1);

Box::Box( double length, double width, double height ) : Hexahedron()
{
//  Make a Box with dimensions length along the x-axis, width along the y-axis, and
//  height along the z-axis.

    constructSurfaces( length, width, height);
}

void Box::resize( double length, double width, double height )
{
//  Resize this Box with dimensions length along the x-axis, width along the y-axis, and
//  height along the z-axis.

    deleteSurfaces();
    constructSurfaces( length, width, height);
}

void Box::constructSurfaces( double length, double width, double height )
{
//  Make 6 planes to describe the box

   max_dimension = 0.5*sqrt( length*length + width*width + height*height);

// Create and add surfaces: all with respect to Box's center
// Note that order above counts here!
    addPlane(  0.5*length*_xhat );
    addPlane( -0.5*length*_xhat );
    addPlane(  0.5*width *_yhat );
    addPlane( -0.5*width *_yhat );
    addPlane(  0.5*height*_zhat );
    addPlane( -0.5*height*_zhat );
}

double Box::length()const { return 2*plane(0).offset();	}
double Box::width() const { return 2*plane(2).offset();}
double Box::height()const { return 2*plane(4).offset();	}

const Vector& Box::xhat()const{ return plane(0).direction();}
const Vector& Box::yhat()const{ return plane(2).direction();}
const Vector& Box::zhat()const{ return plane(4).direction();}

