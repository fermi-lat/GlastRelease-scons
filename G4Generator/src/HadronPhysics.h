#ifndef HadronPhysics_h
#define HadronPhysics_h 1

#include <string>

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4VPhysicsConstructor.hh"
#include "GlastMS/MultipleScatteringFactory.h"

/** 
 * @class HadronPhysics 
 *
 * @brief Hadronic physics processes setup
 *
 *  
 * This class is an implementation of a standard G4VPhysicsConstructor of
 * Geant4. Its main purpouse is to setup a set of particles and physics
 * processes. 
 *
 * In particular, this one activate some hadrons and some relevant hadronic
 * processes
 *  
 * @author F.Longo & F.Paladin
 *    
 * \$Header\$
 */
class HadronPhysics : public G4VPhysicsConstructor
{
  public: 

  //  HadronPhysics(const G4String& name ="hadron");
  
      HadronPhysics(const G4String& name, std::string& physicsChoice, 
          GlastMS::MultipleScatteringFactory& msFactory);
  virtual ~HadronPhysics();
  
 public: 
  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  virtual void ConstructParticle();
  
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();
  
 private:
  
  // Allows to select full hadronic physics
  
  std::string m_physicsChoice;
 
  GlastMS::MultipleScatteringFactory& m_msFactory;
 
};


#endif








