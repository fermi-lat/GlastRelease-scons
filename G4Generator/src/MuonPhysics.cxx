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
#include <iomanip>   

MuonPhysics::MuonPhysics(const G4String& name, 
                         GlastMS::MultipleScatteringFactory& msfactory,
                         GlastMS::EnergyLossFactory& eLossFactory)
                   :  G4VPhysicsConstructor(name), m_msFactory(msfactory), m_eLossFactory(eLossFactory)
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
#include "G4MuNuclearInteraction.hh"


void MuonPhysics::ConstructProcess()
{
  // Purpose and Method: this method is invoked by G4 to build the physics
  //                     processes table

  G4ProcessManager * pManager = 0; 
  
  // Muon Plus Physics

  pManager = G4MuonPlus::MuonPlus()->GetProcessManager();
  
  G4VContinuousDiscreteProcess* fMuPlusIonisation         = 
            m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::IONIZATION);
  G4VContinuousDiscreteProcess* fMuPlusMultipleScattering = m_msFactory(); 
  G4VContinuousDiscreteProcess* fMuPlusBremsstrahlung     = 
            m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::BREMSTRAHLUNG);
  G4MuPairProduction*           fMuPlusPairProduction     = new G4MuPairProduction();
  G4MuNuclearInteraction*       fMuPlusNuclear            = new G4MuNuclearInteraction();

  pManager->AddProcess(fMuPlusMultipleScattering, ordInActive,           1, 1);
  pManager->AddProcess(fMuPlusIonisation,         ordInActive,           2, 2);
  pManager->AddProcess(fMuPlusBremsstrahlung,     ordInActive,           3, 3);
  pManager->AddProcess(fMuPlusPairProduction,     ordInActive,           4, 4);
  pManager->AddProcess(fMuPlusNuclear,            ordInActive, ordInActive, 5);

  // Muon Minus Physics

  pManager = G4MuonMinus::MuonMinus()->GetProcessManager();

  G4VContinuousDiscreteProcess* fMuMinusIonisation         = 
            m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::IONIZATION);
  G4VContinuousDiscreteProcess* fMuMinusMultipleScattering = m_msFactory(); 
  G4VContinuousDiscreteProcess* fMuMinusBremsstrahlung     = 
            m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::BREMSTRAHLUNG);
  G4MuPairProduction*           fMuMinusPairProduction     = new G4MuPairProduction();
  G4MuNuclearInteraction*       fMuMinusNuclear            = new G4MuNuclearInteraction();
  //G4MuonMinusCaptureAtRest* fMuMinusCaptureAtRest = new G4MuonMinusCaptureAtRest();

  pManager->AddProcess(fMuMinusMultipleScattering, ordInActive,           1, 1);
  pManager->AddProcess(fMuMinusIonisation,         ordInActive,           2, 2);
  pManager->AddProcess(fMuMinusBremsstrahlung,     ordInActive,           3, 3);
  pManager->AddProcess(fMuMinusPairProduction,     ordInActive,           4, 4);
  pManager->AddProcess(fMuMinusNuclear,            ordInActive, ordInActive, 5);

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











