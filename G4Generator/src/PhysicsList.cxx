// File and Version Information:
// /$Header/$
//
// Description:
//      PhysicsList .....
//
// Author(s):
//      F.Longo


#include "PhysicsList.h"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ParticleTypes.hh"

#include "GeneralPhysics.h"
#include "EMPhysics.h"
#include "MuonPhysics.h"
#include "HadronPhysics.h"
#include "IonPhysics.h"


PhysicsList::PhysicsList():  G4VModularPhysicsList()
{
  // The default cut value for all particles
  defaultCutValue = 0.1*mm;
  
  // General Physics
  RegisterPhysics( new GeneralPhysics("general") );
  
  // EM Physics
  RegisterPhysics( new EMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics(  new MuonPhysics("muon"));

  // Hadron Physics
  RegisterPhysics(  new HadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics( new IonPhysics("ion"));

}

PhysicsList::~PhysicsList()
{;}

void PhysicsList::SetCuts()
{
  SetCutsWithDefault();   
}






