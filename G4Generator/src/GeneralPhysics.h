#ifndef GeneralPhysics_h
#define GeneralPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

/** 
 * @class GeneralPhysics
 *
 * @brief Generic physics processes setup
 * 
 * This class is an implementation of a standard G4VPhysicsConstructor of
 * Geant4. Its main purpouse is to setup a set of particles and physics
 * processes. 
 *
 * In particular, this one activate the dummy geantino particle and the general
 * decay processes
 *
 *  
 * @author F.Longo 
 *    
 * \$Header\$
 */
class GeneralPhysics : public G4VPhysicsConstructor
{
  public: 
    GeneralPhysics(const G4String& name = "general");
    virtual ~GeneralPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

};


#endif








