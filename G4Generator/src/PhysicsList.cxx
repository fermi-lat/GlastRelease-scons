#include "PhysicsList.h"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ParticleTypes.hh"

PhysicsList::PhysicsList()
{
  defaultCutValue = 0.1*mm;
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
}

PhysicsList::~PhysicsList()
{;}

void PhysicsList::ConstructParticle()
{
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
}

void PhysicsList::ConstructBosons()
{
  G4Gamma::GammaDefinition();
}

void PhysicsList::ConstructLeptons()
{
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

void PhysicsList::ConstructMesons()
{
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
}


void PhysicsList::ConstructBaryons()
{
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}


void PhysicsList::ConstructProcess()
{
  // Define transportation process

  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

void PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
    if (particleName == "gamma") {
      //gamma      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4MultipleScattering(),-1, 1,1);
      pmanager->AddProcess(new G4eIonisation(),       -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung(),   -1,-1,3);      

    } else if (particleName == "e+") {
      //positron      
      pmanager->AddProcess(new G4MultipleScattering(),-1, 1,1);
      pmanager->AddProcess(new G4eIonisation(),       -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung(),   -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation(),  0,-1,4);      

    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
     //muon  
     pmanager->AddProcess(new G4MultipleScattering(),-1, 1,1);
     pmanager->AddProcess(new G4MuIonisation(),      -1, 2,2);
     pmanager->AddProcess(new G4MuBremsstrahlung(),  -1,-1,3);
     pmanager->AddProcess(new G4MuPairProduction(),  -1,-1,4);       
     
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4hIonisation(),       -1,2,2);      
    } 
  }
}

#include "G4Decay.hh"

void PhysicsList::ConstructGeneral()
{
  // Add Decay Process
   G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

void PhysicsList::SetCuts()
{
  SetCutsWithDefault();   
}


