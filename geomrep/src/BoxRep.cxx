// $Header$
//  Author: Toby Burnett
//
// Display of any Hexahedron

#include "geomrep/BoxRep.h"

#include "geometry/Hexahedron.h"

BoxRep::BoxRep(const Shape& shape)
    : m_box( dynamic_cast<const Hexahedron& >(shape) )
{};
    

void BoxRep::update()
{
    const int lines[]=
    {	0,2,3,1,0, 
        4,6,7,5,4, -1,
        2,6,-1,
	3,7,-1,
        1,5,-1
    };

    for(int i=0; i<8; i++)  {
	addVertex(m_box.vertex(i));
    }

    indexedLineSet(lines, sizeof(lines)/sizeof(int));
}


