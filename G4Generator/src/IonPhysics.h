#ifndef IonPhysics_h
#define IonPhysics_h 1

#include <string>

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "Geant4/MultipleScatteringFactory.h"

/** 
 * @class IonPhysics 
 *
 * @brief Ions physics processes setup
 *
 * This class is an implementation of a standard G4VPhysicsConstructor of
 * Geant4. Its main purpouse is to setup a set of particles and physics
 * processes. 
 *
 * In particular, this one activate the Geant4 ions and their processes
 *  
 * @author F.Longo & F.Paladin
 *    
 * \$Header\$
 */

class IonPhysics : public G4VPhysicsConstructor
{
 public: 
  
  IonPhysics(const G4String& name, std::string& physicsChoice, 
      Geant4::MultipleScatteringFactory& msfactory);
  virtual ~IonPhysics();
  
 public: 
  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  virtual void ConstructParticle();
  
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();
  
 private:
  
  // Allow to select full hadronic physics
  
  std::string m_physicsChoice;
    Geant4::MultipleScatteringFactory& m_msFactory;

};


#endif

