// $Id$
//
//
#ifndef __HEXAHEDRON_H
#define __HEXAHEDRON_H

#include "geometry/Volume.h"
#include "geometry/Plane.h"

/**
// Define the abstract class Hexahedron, a solid bounded by six plane suraces.
<PRE>
Vertex numbering:
          			 ^ y
          			 |
		   	1/5-------------0/4
			 |	 |       |
			 |	 |       |
	       	 --------|-------+-------|------> x
			 |	 |       |
			 |	 |       |
			3/7------+------2/6
				 |
     numbers are indeces for +z/-z vertices
</PRE>
*/
class Hexahedron: public Volume
{
public:
    
    Hexahedron();
    
    
    /// return the vertex number n	(see fig for indexing)
    Point vertex(int n)const;
    
protected:
    
    /// simplify construction
    void addPlane(const Vector& v){addSurface(new Plane(center(),v));}
    
    /// avoid having to downcast in sub classes
    const Plane& plane(int i)const{return (const Plane&)Volume::surface(i);}
};

#endif

