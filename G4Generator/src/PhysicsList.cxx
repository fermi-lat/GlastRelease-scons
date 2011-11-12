// File and Version Information:
// $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/G4Generator/src/PhysicsList.cxx,v 1.24 2011/08/26 23:59:15 flongo Exp 
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

// GLAST

// others
#include "EmStandardPhysics.hh"
#include "EmLowEnergyPhysics.hh"

#include "EpaxIonPhysics.h"

//#include "HadronPhysicsLHEP.hh"
//#include "HadronPhysicsQGSP.hh"
//#include "HadronPhysicsQGSP_BIC.hh"
#include "HadronPhysicsQGSP_BERT.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4DecayPhysics.hh"
//#include "G4IonPhysics.hh"


// new 

#include "G4QStoppingPhysics.hh" 
#include "G4HadronElasticPhysicsLHEP.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4IonBinaryCascadePhysics.hh"



#include <cstdlib>

PhysicsList::PhysicsList(double cutValue, const std::string& physicsChoice, 
                         const std::string& physicsTable, const std::string& physicsDir)
                         :  G4VModularPhysicsList()
{
  // The default cut value for all particles
  defaultCutValue = cutValue;

  // Physics Choice

  m_physicsChoice = physicsChoice;

  // Physics Tables

  m_physicsTable = physicsTable;
  m_physicsDir = physicsDir;


 std::cout << physicsChoice << std::endl; // this is for hadronic stuff

  //          EmLow
  //          EmStd
  //          EmLiv // New
  //          QGSP_BERT_EPAX // New

 if (m_physicsChoice=="EmLow")
   {
     G4int ver=0;
  
     // EM Physics
     
     RegisterPhysics( new EmLowEnergyPhysics("lowEnergy EM", ver));
     
     // Decay Physics - i.e. decay
     
     RegisterPhysics( new G4DecayPhysics("decay") );
     
     // Ion Physics
     
     RegisterPhysics( new EpaxIonPhysics("ion", m_physicsChoice));
     
     RegisterPhysics(  new HadronPhysicsQGSP_BERT("hadron")); // to be updated
 
   }
 if (m_physicsChoice=="EmLiv")
   {
     G4int ver=0;
  
     // EM Physics
     
     RegisterPhysics( new G4EmLivermorePhysics(ver));
     RegisterPhysics( new G4EmExtraPhysics("extra EM"));
     
     // Decay Physics - i.e. decay
     
     RegisterPhysics( new G4DecayPhysics("decay") );
     
     // Ion Physics
     
     RegisterPhysics( new EpaxIonPhysics("ion", m_physicsChoice));

     // Hadron Phyiscs
     RegisterPhysics( new G4NeutronTrackingCut(ver) );
     RegisterPhysics(  new HadronPhysicsQGSP_BERT("hadron")); // to be updated
     RegisterPhysics(  new G4HadronElasticPhysics("hadronElast"));
     RegisterPhysics( new G4QStoppingPhysics("stopping"));
     
   }
 if (m_physicsChoice=="EmStd") 
   {

     G4int ver=0;
  
     // EM Physics
     
     RegisterPhysics( new EmStandardPhysics("standard EM", ver));

     // Synchroton Radiation & GN Physics
  
     RegisterPhysics( new G4EmExtraPhysics("extra EM"));

     // Decay Physics - i.e. decay
     
     RegisterPhysics( new G4DecayPhysics("decay") );
     
     // Full Ion Physics
     
     RegisterPhysics( new EpaxIonPhysics("ion", m_physicsChoice));


     // Neutrons 

     RegisterPhysics( new G4NeutronTrackingCut(ver) );
     
     RegisterPhysics(  new HadronPhysicsQGSP_BERT("hadron"));
     
     // Hadron Elastic Phys
     
     RegisterPhysics(  new G4HadronElasticPhysics("hadronElast"));
     
     // Stopping Physics
     
     RegisterPhysics( new G4QStoppingPhysics("stopping"));
     
     
   }
 if (m_physicsChoice=="QGSP_BERT_EPAX") 
   {
     G4int ver=0;
     
     // EM standard physics (default)
     
     RegisterPhysics( new G4EmStandardPhysics(ver));
     
     // Synchroton Radiation & GN Physics
     
     RegisterPhysics( new G4EmExtraPhysics("extra EM"));

     // Decay Physics - i.e. decay
     
     RegisterPhysics( new G4DecayPhysics("decay") );
     
     // Ion Physics
    
     RegisterPhysics( new EpaxIonPhysics("ion", m_physicsChoice));
     
     // Neutrons 

     RegisterPhysics( new G4NeutronTrackingCut(ver) );
     
     // Hadron Physics
     
     RegisterPhysics(  new HadronPhysicsQGSP_BERT("hadron"));
     
     // Hadron Elastic Phys
     
     RegisterPhysics(  new G4HadronElasticPhysics("hadronElast"));
     
     // Stopping Physics
     
     RegisterPhysics( new G4QStoppingPhysics("stopping"));
     
     
   }




}
PhysicsList::~PhysicsList()
{;}


void PhysicsList::SetCuts()
{
#ifndef SCons
  std::string physics_basedir = getenv("G4GENERATORROOT");
#else
  // Shouldn't matter since this routine is normally called with
  // m_physicsTable = "build"
  std::string physics_basedir("stuffAndNonsense");
#endif

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
