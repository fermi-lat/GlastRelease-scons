//      $Id$
//   Author: T. Burnett
//  Project: Arve graphics
//
// Representation of a Cone with ArveScene

#ifndef CONEREP_H
#define CONEREP_H
#include "gui/DisplayRep.h"

class Cone;

class ConeRep : public gui::DisplayRep {

public:
    ConeRep(const Cone& cone):m_cone(cone){};
    void update();
private:
    const Cone& m_cone;
};
#endif //CONEREP_H
