// $Header$
//

#ifndef __INTERSECTION_H
#define __INTERSECTION_H

#include "geometry/Ray.h"
#include "geometry/Surface.h"

/**
Special class that handles the association of a Ray and a
Surface object to find the intersection if any.
<BR>
Both Ray and Surface are polymorphic
*/
class Intersection 
{
public:
    
    Intersection(const Ray& r, const Surface & s )
        : ray(r), surf(s)
    {}
    
    
    /// return the actual intersection, or FLT_MAX. sign>0 -> leaving; sign<0  -> entering
    double distance(double dist, int sign);
    
private:
    const Ray& ray;
    const Surface& surf;
    
};
#endif

