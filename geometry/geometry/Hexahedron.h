// $Id$
//
//
// Define the abstract class Hexahedron, a solid bounded by six plane
// surfaces.  Subclasses are Box, Trap, Trd1, Trd2 
//
//          			 ^ y
//          			 |
//		   	1/5-------------0/4
//			 |	 |       |
//			 |	 |       |
//	       	 --------|-------+-------|------> x
//			 |	 |       |
//			 |	 |       |
//			3/7------+------2/6
//				 |
//     numbers are indeces for +z/-z vertices
#ifndef __HEXAHEDRON_H
#define __HEXAHEDRON_H

#include "geometry/Volume.h"
#include "geometry/Plane.h"

class Hexahedron: public Volume
{
public:
    
    Hexahedron();
    
    
    
    Point vertex(int n)const;
    // return the vertex number n	(see fig for indexing)
    
protected:
    
    void addPlane(const Vector& v){addSurface(new Plane(center(),v));}
    // simplify construction
    
    const Plane& plane(int i)const{return (const Plane&)Volume::surface(i);}
    // avoid having to downcast in sub classes
};

#endif

