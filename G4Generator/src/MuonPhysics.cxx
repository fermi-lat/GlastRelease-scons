// File and Version Information:
// $Header$
//
// Description: This class manages the building of muon/tau and their processes
//
// Author(s):
//      F.Longo

#include "MuonPhysics.h"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   

MuonPhysics::MuonPhysics(const G4String& name, 
                         GlastMS::MultipleScatteringFactory& msfactory)
                   :  G4VPhysicsConstructor(name), m_msFactory(msfactory)
{
}

MuonPhysics::~MuonPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"

void MuonPhysics::ConstructParticle()
{
  // Purpose and Method: this method is invoked by G4 to build the particles
  //                     classes

  // Mu

  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  // Tau
  G4TauMinus::TauMinusDefinition();
  G4TauPlus::TauPlusDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();

}


#include "G4ProcessManager.hh"

//#include "G4MultipleScattering.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuIonisation.hh"
#include "G4hIonisation.hh"
#include "G4MuonMinusCaptureAtRest.hh"


void MuonPhysics::ConstructProcess()
{
  // Purpose and Method: this method is invoked by G4 to build the physics
  //                     processes table

  G4ProcessManager * pManager = 0; 
  
  // Muon Plus Physics

  pManager = G4MuonPlus::MuonPlus()->GetProcessManager();
  
  G4VContinuousDiscreteProcess*   fMuPlusMultipleScattering = m_msFactory(); 
  G4MuBremsstrahlung*     fMuPlusBremsstrahlung = new G4MuBremsstrahlung();
  G4MuPairProduction*     fMuPlusPairProduction = new G4MuPairProduction();
  G4MuIonisation*         fMuPlusIonisation = new G4MuIonisation();
  pManager->AddProcess(fMuPlusIonisation, ordInActive,2, 2);
  
  pManager->AddDiscreteProcess(fMuPlusBremsstrahlung);
  pManager->AddDiscreteProcess(fMuPlusPairProduction);
  pManager->AddProcess(fMuPlusMultipleScattering);
  pManager->SetProcessOrdering(fMuPlusMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fMuPlusMultipleScattering, idxPostStep,  1);

  // Muon Minus Physics

  pManager = G4MuonMinus::MuonMinus()->GetProcessManager();

  G4VContinuousDiscreteProcess*   fMuMinusMultipleScattering = m_msFactory();
  G4MuBremsstrahlung*     fMuMinusBremsstrahlung = new G4MuBremsstrahlung();
  G4MuPairProduction*     fMuMinusPairProduction = new G4MuPairProduction();
  G4MuIonisation*         fMuMinusIonisation = new  G4MuIonisation();
  //G4MuonMinusCaptureAtRest* fMuMinusCaptureAtRest = new G4MuonMinusCaptureAtRest();
  pManager->AddProcess(fMuMinusIonisation, ordInActive,2, 2);
  pManager->AddDiscreteProcess(fMuMinusBremsstrahlung);
  pManager->AddDiscreteProcess(fMuMinusPairProduction);
  pManager->AddProcess(fMuMinusMultipleScattering);
  pManager->SetProcessOrdering(fMuMinusMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fMuMinusMultipleScattering, idxPostStep,  1);
  //  pManager->AddRestProcess(fMuMinusCaptureAtRest);

  // Tau Plus Physics

  pManager = G4TauPlus::TauPlus()->GetProcessManager();

  G4VContinuousDiscreteProcess*   fTauPlusMultipleScattering = m_msFactory();
  G4hIonisation*          fTauPlusIonisation = new G4hIonisation();
  pManager->AddProcess(fTauPlusIonisation, ordInActive,2, 2);
  pManager->AddProcess(fTauPlusMultipleScattering);
  pManager->SetProcessOrdering(fTauPlusMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fTauPlusMultipleScattering, idxPostStep,  1);

  // Tau Minus Physics

  pManager = G4TauMinus::TauMinus()->GetProcessManager();

  G4VContinuousDiscreteProcess*   fTauMinusMultipleScattering =  m_msFactory();
  G4hIonisation*          fTauMinusIonisation = new G4hIonisation();
  pManager->AddProcess(fTauMinusIonisation, ordInActive,2, 2);
  pManager->AddProcess(fTauMinusMultipleScattering);
  pManager->SetProcessOrdering(fTauMinusMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fTauMinusMultipleScattering, idxPostStep,  1);

}











