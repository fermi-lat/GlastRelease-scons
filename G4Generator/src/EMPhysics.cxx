// File and Version Information:
// $Header$
//
// Description: This class manages the building of gamma/electron/positron and
// their processes
//
// Author(s):
//      F.Longo

#include "EMPhysics.h"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   


EMPhysics::EMPhysics(const G4String& name, GlastMS::MultipleScatteringFactory& msFactory)
               :  G4VPhysicsConstructor(name), m_msFactory(msFactory)
{
}

EMPhysics::~EMPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

void EMPhysics::ConstructParticle()
{
  // Purpose and Method: this method is invoked by G4 to build the particles
  //                     classes

  // gamma
  G4Gamma::GammaDefinition();
 
  // electron
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
}


#include "G4ProcessManager.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
//THB #include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

void EMPhysics::ConstructProcess()
{
  // Purpose and Method: this method is invoked by G4 to build the physics
  //                     processes table

  G4ProcessManager * pManager = 0;

  // Gamma Physics

  pManager = G4Gamma::Gamma()->GetProcessManager();
  
  G4PhotoElectricEffect* thePhotoEffect = new G4PhotoElectricEffect();
  G4ComptonScattering* theComptonEffect = new G4ComptonScattering();
  G4GammaConversion* thePairProduction = new  G4GammaConversion();
  pManager->AddDiscreteProcess(thePhotoEffect);
  pManager->AddDiscreteProcess(theComptonEffect);
  pManager->AddDiscreteProcess(thePairProduction);

  // Electron Physics

  pManager = G4Electron::Electron()->GetProcessManager();

  G4VContinuousDiscreteProcess* theElectronMultipleScattering = m_msFactory();
  G4eIonisation* theElectronIonisation = new G4eIonisation();
  G4eBremsstrahlung* theElectronBremsStrahlung = new G4eBremsstrahlung();
  pManager->AddDiscreteProcess(theElectronBremsStrahlung);  
  pManager->AddProcess(theElectronIonisation, ordInActive,2, 2);
  pManager->AddProcess(theElectronMultipleScattering);
  pManager->SetProcessOrdering(theElectronMultipleScattering, 
                               idxAlongStep,  1);
  pManager->SetProcessOrdering(theElectronMultipleScattering, 
                               idxPostStep,  1);

  //Positron Physics

  pManager = G4Positron::Positron()->GetProcessManager();

  G4VContinuousDiscreteProcess* thePositronMultipleScattering  =  m_msFactory();
  G4eIonisation* thePositronIonisation = new  G4eIonisation(); 
  G4eBremsstrahlung* thePositronBremsStrahlung = new G4eBremsstrahlung();  
  G4eplusAnnihilation* theAnnihilation = new G4eplusAnnihilation();
  pManager->AddDiscreteProcess(thePositronBremsStrahlung);
  pManager->AddDiscreteProcess(theAnnihilation);
  pManager->AddRestProcess(theAnnihilation);
  pManager->AddProcess(thePositronIonisation, ordInActive,2, 2);
  pManager->AddProcess(thePositronMultipleScattering);
  pManager->SetProcessOrdering(thePositronMultipleScattering, 
                               idxAlongStep,  1);
  pManager->SetProcessOrdering(thePositronMultipleScattering, 
                               idxPostStep,  1);
}



