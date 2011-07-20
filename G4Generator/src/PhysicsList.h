#ifndef PhysicsList_h
#define PhysicsList_h 1

#include <string>

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

/** 
 * @class PhysicsList
 *
 * @brief Physics list class for G4Generator
 *
 * This class defines and activates all the physics processes used by
 * G4Generator: in particular it sets
 * 
 * - the particles to be used in the simulation
 * - the range cuts for each particle
 * - the physics processes to be simulated
 *  
 * @author F.Longo 
 *    
 * $Header$
 */
class PhysicsList: public G4VModularPhysicsList
{
 public:
  PhysicsList(double cutValue, const std::string& physicsChoice, 
      const std::string& physicsTable, const std::string& physicsDir);
  ~PhysicsList();
  
 public:

  // Construct particle and physics process
  // virtual void ConstructParticle();
  // virtual void ConstructProcess();
  
  /// This method set all the physics cuts for the simulation
  virtual void SetCuts();
  

 private:
  
  std::string m_physicsChoice;

  std::string m_physicsTable; // still needed? 
  std::string m_physicsDir;

  /*  
      G4VPhysicsConstructor* m_GeneralPhysics;
      G4VPhysicsConstructor* m_EMPhysics;
      G4VPhysicsConstructor* m_MuonPhysics;
      G4VPhysicsConstructor* m_HadronPhysics;
      G4VPhysicsConstructor* m_IonPhysics;
  */

};

#endif







