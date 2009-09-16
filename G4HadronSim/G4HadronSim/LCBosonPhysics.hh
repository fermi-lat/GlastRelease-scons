//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Title:  Boson physics for a Linear Collider Detector                   //
//  Date:   16 June 2004                                                   //
//  Author: D.H. Wright (SLAC)                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

#ifndef LCBosonPhysics_h
#define LCBosonPhysics_h 1

#include "G4VPhysicsConstructor.hh"

class LCBosonPhysics : public G4VPhysicsConstructor
{
  public: 
    LCBosonPhysics(const G4String& name = "boson");
    virtual ~LCBosonPhysics();

    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

};

#endif








