// File and Version Information:
// $Header$
//
// Description: This class manages the building of particles definitions and
// physics processes setup by creating a set of specialized classes and
// registering them
//      
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

#include <cstdlib>

PhysicsList::PhysicsList(double cutValue, const std::string& physicsChoice, 
                         const std::string& physicsTable, const std::string& physicsDir,
                         GlastMS::MultipleScatteringFactory& msFactory
                         ):  G4VModularPhysicsList()
{
  // The default cut value for all particles
  defaultCutValue = cutValue;

  // Physics Choice

  m_physicsChoice = physicsChoice;

  // Physics Tables

  m_physicsTable = physicsTable;
  m_physicsDir = physicsDir;

  
  // General Physics
  RegisterPhysics( new GeneralPhysics("general") );
  
  // EM Physics 

  RegisterPhysics( new EMPhysics("standard EM", msFactory));

  // Muon Physics

  RegisterPhysics(  new MuonPhysics("muon", msFactory));

  // Full or EM Hadron Physics

  RegisterPhysics(  new HadronPhysics("hadron", m_physicsChoice, msFactory));
  // RegisterPhysics(  new HadronPhysics("hadron"));
  
  // Full or EM Ion Physics

  RegisterPhysics( new IonPhysics("ion", m_physicsChoice, msFactory));
  //RegisterPhysics( new IonPhysics("ion"));
  

}

PhysicsList::~PhysicsList()
{;}

void PhysicsList::SetCuts()
{

  std::string physics_basedir = getenv("G4GENERATORROOT");
  std::string physics_Dir = physics_basedir + "/" + m_physicsDir;
  if (m_physicsTable=="retrieve")
    {  
      SetPhysicsTableRetrieved(physics_Dir);
    }
  SetCutsWithDefault(); 
  if (m_physicsTable=="store")
    {  
      if(StorePhysicsTable(physics_Dir))
	{// std::cout << "done Physics Store" << std::endl;
	}
    }
}














