// $Id$
// This file is part of Gismo 2
//
//
// Intersection is a special class that handles the association of a Ray and a
// Surface object to find the intersection if any.
//
// Both Ray and Surface are polymorphic
//
// Actual numerical solution is performed by the superclass RootFinder


#ifndef __INTERSECTION_H
#define __INTERSECTION_H

#include "geometry/Ray.h"
#include "geometry/Surface.h"


class Intersection 
{
public:

    Intersection(const Ray& r, const Surface & s )
    : ray(r), surf(s)
    {}

  
    double distance(double dist, int sign);
    // return the actual intersection, or FLT_MAX. sign>0 -> leaving; sign<0  -> entering

private:
    const Ray& ray;
    const Surface& surf;

};
#endif

