#ifndef EMPhysics_h
#define EMPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "GlastMS/MultipleScatteringFactory.h"


/** 
 * @class EMPhysics 
 *
 * @brief Electromagnetic physics processes setup
 *
 * This class is an implementation of a standard G4VPhysicsConstructor of
 * Geant4. Its main purpouse is to setup a set of particles and physics
 * processes. 
 *
 * In particular, this one activate the gamma, the electron, the positron and
 * the electromagnetic processes related to them
 *  
 * @author F.Longo 
 *    
 * \$Header\$
 */
class EMPhysics : public G4VPhysicsConstructor
{
  public: 
      EMPhysics(const G4String& name , GlastMS::MultipleScatteringFactory& msfactory);
    virtual ~EMPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

private:
    GlastMS::MultipleScatteringFactory& m_msFactory;
};


#endif





