//     $Id$
// Project: Atlas Reconstruction prototype
//  Author: Toby Burnett
//
//  Base class for atlas digitizer (or G4's G4SensitiveDetector )
//
// Inherits from CoordSystem to provide for an internal geometry

#ifndef DIGITIZER_H
#define DIGITIZER_H

#include "geometry/CoordSystem.h"
#include "geometry/Ray.h"

class Digitizer : public CoordSystem {
public:

    virtual void addHit(const Ray& pt, float eloss)=0;
    // basic interface to the simulation: receive a hit in the form
    // of energy loss from a charged particle at the given point

    virtual void clear() {}
    virtual void generateResponse() {}
    // these, if implemented, handle accumulation of data

    virtual ~Digitizer(){}
protected:
    Digitizer(){};
    // abstract class
};
#endif
