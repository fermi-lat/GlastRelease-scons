#ifndef MuonPhysics_h
#define MuonPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "Geant4/MultipleScatteringFactory.h"

/** 
 * @class MuonPhysics 
 *
 * @brief Muons physics processes setup
 *
 * This class is an implementation of a standard G4VPhysicsConstructor of
 * Geant4. Its main purpouse is to setup a set of particles and physics
 * processes. 
 *
 * In particular, this one activate the muon and tau particles and their
 * processes
 *  
 * @author F.Longo 
 *    
 * \$Header\$
 */
class MuonPhysics : public G4VPhysicsConstructor
{
  public: 
    MuonPhysics(const G4String& name, 
        Geant4::MultipleScatteringFactory& msfactory);
    virtual ~MuonPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();
private:
    Geant4::MultipleScatteringFactory& m_msFactory;

};


#endif

