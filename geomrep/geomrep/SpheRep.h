// $Header$
//  Author: T. Burnett
// Project: Arve graphics
//
// Representation of a Tube with ArveScene

#ifndef TUBEREP_H
#define TUBEREP_H

#include "gui/DisplayRep.h"

class Sphe;

class SpheRep : public gui::DisplayRep {

public:
    SpheRep(const Sphe& s):m_sphe(s){};
    void update();
    const Sphe& returnSphe(){return m_sphe;}
private:
    const Sphe& m_sphe;
};
#endif //SPHEREP_H

