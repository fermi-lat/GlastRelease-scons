// $Header$
//  Author: T. Burnett
// Project: Arve graphics
//
// Representation of a Tube with ArveScene

#ifndef TUBEREP_H
#define TUBEREP_H

#include "graphics/ArvePlottableRep.h"
#include "geometry/Tube.h"

class Shape;

class TubeRep : public ArvePlottableRep {
    
public:
    TubeRep(const Tube& tube):m_tube(tube){};
    
    TubeRep(const Shape& shape); 
    // special down-casting constructor
    
    void update();
private:
    const Tube& m_tube;
};
#endif //TUBEREP_H

