// $Header$
//  Author: T. Burnett
// Project: Arve graphics
//
// Representation of a Tube with ArveScene

#ifndef TUBEREP_H
#define TUBEREP_H

#include "graphics/ArvePlottableRep.h"

class Sphe;

class SpheRep : public ArvePlottableRep {

public:
    SpheRep(const Sphe& s):m_sphe(s){};
    void update();
private:
    const Sphe& m_sphe;
};
#endif //SPHEREP_H

