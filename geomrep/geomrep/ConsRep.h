// $Heading: ConeRep.h $
//   Author: T. Burnett
//
// Representation of a Cons with ArveScene

#ifndef CONSREP_H
#define CONSREP_H
#include "gui/DisplayRep.h"

class Cons;

class ConsRep : public gui::DisplayRep {

public:
    ConsRep(const Cons& cons):m_cons(cons){};
    void update();
private:
    const Cons& m_cons;
};
#endif //CONSREP_H

