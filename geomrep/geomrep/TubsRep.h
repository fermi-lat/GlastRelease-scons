// $Header$
//  Author: T. Burnett
// Project: Arve graphics
//
// Representation of a Tube with ArveScene

#ifndef TUBSREP_H
#define TUBSREP_H
#include "graphics/ArvePlottableRep.h"

class Tubs;

class TubsRep : public ArvePlottableRep {

public:
    TubsRep(const Tubs& tubs):m_tubs(tubs){};
    void update();
private:
    const Tubs& m_tubs;
};
#endif //TUBSREP_H

