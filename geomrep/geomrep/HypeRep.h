// $Header$
//  Author: T. Burnett
// Project: Arve graphics
//
// Representation of a Hype with ArveScene

#ifndef TUBSREP_H
#define HYPEREP_H
#include "graphics/ArvePlottableRep.h"

class Hype;

class HypeRep : public ArvePlottableRep {

public:
    HypeRep(const Hype& hype):m_hype(hype){};
    void update();
private:
    const Hype& m_hype;
};
#endif //HYPEREP_H

