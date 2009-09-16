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

// GLAST

#include "GeneralPhysics.h"
#include "EMPhysics.h"
#include "MuonPhysics.h"
#include "HadronPhysics.h"
#include "IonPhysics.h"

// LC 

#include "G4HadronSim/LCBosonPhysics.hh"
#include "G4HadronSim/LCDecayPhysics.hh"
#include "G4HadronSim/LCHadronPhysics.hh"
#include "G4HadronSim/LCIonPhysics.hh"
#include "G4HadronSim/LCLeptonPhysics.hh"

// Space

#include "G4HadronSim/SEBosonPhysics.hh"
#include "G4HadronSim/SEDecayPhysics.hh"
#include "G4HadronSim/SEHadronPhysics.hh"
#include "G4HadronSim/SEIonPhysics.hh"
#include "G4HadronSim/SELeptonPhysics.hh"
#include "G4HadronSim/SENeutronPhysics.hh"

// others

#include "G4HadronSim/HadronPhysicsLHEP.hh"
#include "G4HadronSim/HadronPhysicsLHEP_BIC.hh"
#include "G4HadronSim/HadronPhysicsLHEP_BERT.hh"
#include "G4HadronSim/HadronPhysicsQGSP.hh"
#include "G4HadronSim/HadronPhysicsQGSP_BIC.hh"
#include "G4HadronSim/HadronPhysicsQGSP_BERT.hh"
#include "G4HadronSim/HadronPhysicsQGSC.hh"
#include "G4HadronSim/HadronPhysicsQGSC_LEAD.hh"
#include "G4HadronSim/G4EmStandardPhysics.hh"
#include "G4HadronSim/G4EmExtraPhysics.hh"
#include "G4HadronSim/G4EmLowEnergyPhysics.hh"
#include "G4HadronSim/G4DecayPhysics.hh"
#include "G4HadronSim/G4IonPhysics.hh"
#include "EpaxIonPhysics.h"


#include <cstdlib>

PhysicsList::PhysicsList(double cutValue, const std::string& physicsChoice, 
                         const std::string& physicsTable, const std::string& physicsDir,
                         GlastMS::MultipleScatteringFactory& msFactory,
                         GlastMS::EnergyLossFactory& eLossFactory
                         ):  G4VModularPhysicsList()
{
  // The default cut value for all particles
  defaultCutValue = cutValue;

  // Physics Choice

  m_physicsChoice = physicsChoice;

  // Physics Tables

  m_physicsTable = physicsTable;
  m_physicsDir = physicsDir;


 std::cout << physicsChoice << std::endl; // this is for hadronic stuff

  // options: GLAST
  //          LHEP
  //          LHEP_BIC
  //          LHEP_BERT
  //          QGSP
  //          QGSP_BIC
  //          QGSP_BERT
 //           QGSP_BERT_LPM
  //          LC 
  //          Space
 //           LowEnergy 

 if (m_physicsChoice=="GLAST")
   {
     // General Physics
     RegisterPhysics( new GeneralPhysics("general") );
     
     // EM Physics 
     
     RegisterPhysics( new EMPhysics("standard EM", msFactory, eLossFactory));
     
     // Muon Physics
     
     RegisterPhysics(  new MuonPhysics("muon", msFactory, eLossFactory));

     // Full Hadron Physics
     
     RegisterPhysics(  new HadronPhysics("hadron", m_physicsChoice, msFactory));
     // RegisterPhysics(  new HadronPhysics("hadron"));
     
     // Full Ion Physics
     
     RegisterPhysics( new IonPhysics("ion", m_physicsChoice, msFactory));
     
     //     // Full Ion Physics
     
     //RegisterPhysics( new IonPhysics("ion", msFactory));
   }
 if (m_physicsChoice=="LC")
   {
     
     RegisterPhysics( new LCDecayPhysics("decay"));
     
     // Bosons (gamma + geantinos)
     RegisterPhysics( new LCBosonPhysics("boson"));
      
     // Leptons
     RegisterPhysics( new LCLeptonPhysics("lepton",msFactory, eLossFactory));
     
     // Hadron Physics
     RegisterPhysics( new LCHadronPhysics("hadron", msFactory));
     
      // Ion Physics
     RegisterPhysics( new LCIonPhysics("ion", msFactory));
   }
 if (m_physicsChoice=="Space")
   {
     RegisterPhysics( new SEDecayPhysics("decay"));
     
     // Bosons (gamma + geantinos)
     RegisterPhysics( new SEBosonPhysics("boson"));
      
     // Leptons
     RegisterPhysics( new SELeptonPhysics("lepton",msFactory, eLossFactory));
     
     // Hadron Physics
     RegisterPhysics( new SEHadronPhysics("hadron",msFactory));

     // Neutron Physics
     RegisterPhysics( new SENeutronPhysics("neutron"));
     
      // Ion Physics
     RegisterPhysics( new SEIonPhysics("ion",msFactory));
   }
 if (m_physicsChoice=="LowEnergy")
   {
     G4int ver=1;
  
     // EM Physics
     
     RegisterPhysics( new G4EmLowEnergyPhysics("lowEnergy EM", ver, msFactory, eLossFactory));
     
     // Decay Physics - i.e. decay
     
     RegisterPhysics( new G4DecayPhysics("decay",ver) );
     
     // Ion Physics
     
     RegisterPhysics( new G4IonPhysics("ion")); // to be updated

     
     RegisterPhysics(  new HadronPhysicsLHEP_BERT("hadron")); // to be updated
 
   }
 if (m_physicsChoice!="GLAST"&&m_physicsChoice!="LC"&&m_physicsChoice!="Space"&&m_physicsChoice!="LowEnergy") 
   {

     G4int ver=1;
     if (m_physicsChoice=="QGSP_BERT_LPM") ver=2;
  
     // EM Physics
     
     RegisterPhysics( new G4EmStandardPhysics("standard EM", ver, msFactory, eLossFactory));
     
     // Synchroton Radiation & GN Physics
  
     RegisterPhysics( new G4EmExtraPhysics("extra EM"));

     // Decay Physics - i.e. decay
     
     RegisterPhysics( new G4DecayPhysics("decay",ver) );
     
     // Ion Physics
     //RegisterPhysics( new G4IonPhysics("ion"));

     // Full Ion Physics
     
     RegisterPhysics( new EpaxIonPhysics("ion", m_physicsChoice, msFactory));


     if (m_physicsChoice=="LHEP")
     // Hadron Physics
       RegisterPhysics(  new HadronPhysicsLHEP("hadron"));

     if (m_physicsChoice=="LHEP_BIC")
       RegisterPhysics(  new HadronPhysicsLHEP_BIC("hadron"));

     if (m_physicsChoice=="LHEP_BERT")
       RegisterPhysics(  new HadronPhysicsLHEP_BERT("hadron"));

     if (m_physicsChoice=="QGSP")
       RegisterPhysics(  new HadronPhysicsQGSP("hadron"));

     if (m_physicsChoice=="QGSP_BERT")
       RegisterPhysics(  new HadronPhysicsQGSP_BERT("hadron"));

     if (m_physicsChoice=="QGSP_BERT_LPM")
       RegisterPhysics(  new HadronPhysicsQGSP_BERT("hadron"));

     if (m_physicsChoice=="QGSP_BIC")
       RegisterPhysics(  new HadronPhysicsQGSP_BIC("hadron"));

     if ( m_physicsChoice == "QGSC" )
         RegisterPhysics(new HadronPhysicsQGSC("hadron"));

     if ( m_physicsChoice == "QGSC_LEAD" )
         RegisterPhysics(new HadronPhysicsQGSC_LEAD("hadron"));

   }
}
PhysicsList::~PhysicsList()
{;}

/*void PhysicsList::ConstructParticle()
{
  m_GeneralPhysics->ConstructParticle();
  m_EMPhysics->ConstructParticle();
  m_MuonPhysics->ConstructParticle();
  m_HadronPhysics->ConstructParticle();
  m_IonPhysics->ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
  AddTransportation();

  m_GeneralPhysics->ConstructProcess();
  m_EMPhysics->ConstructProcess();
  m_MuonPhysics->ConstructProcess();
  m_HadronPhysics->ConstructProcess();
  m_IonPhysics->ConstructProcess();
}
*/


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
