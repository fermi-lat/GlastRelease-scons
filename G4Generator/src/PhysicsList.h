#ifndef PhysicsList_h
#define PhysicsList_h 1

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
 * \$Header\$
 */
class PhysicsList: public G4VModularPhysicsList
{
 public:
  PhysicsList();
  ~PhysicsList();
  
 public:
  
  /// This method set all the physics cuts for the simulation
  virtual void SetCuts();
};

#endif







