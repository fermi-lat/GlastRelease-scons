//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Title:  Neutron Physics for Space Electronics Applications             //
//  Date:   20 March 2005                                                  //
//  Author: D.H. Wright (SLAC)                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

#ifndef SENeutronPhysics_h
#define SENeutronPhysics_h 1

#include "G4VPhysicsConstructor.hh"

class SENeutronPhysics : public G4VPhysicsConstructor
{
  public: 
    SENeutronPhysics(const G4String& name ="neutron");
    virtual ~SENeutronPhysics();

    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

};

#endif





