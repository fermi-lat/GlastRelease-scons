//  $Header$
//   Author: Toby Burnett
//  Project: Arve graphics 
// 
// Display of any Hexahedron by ArveScene

#ifndef BOXREP_H
#define BOXREP_H

#include "graphics/ArvePlottableRep.h"
class Hexahedron;
class Shape;

class BoxRep : public ArvePlottableRep {
    // Represent a Hexahedron with Arve graphics
public:
    BoxRep(const Hexahedron& box) 
	: m_box(box){};

    BoxRep(const Shape& shape); 
    // special down-casting constructor
    
    virtual void update();

private:
    const Hexahedron& m_box;

};
#endif //BOXREP_H
