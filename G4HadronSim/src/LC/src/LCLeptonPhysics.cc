//
////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Title:  Lepton physics for a Linear Collider Detector                 //
//  Date:   17 June 2004                                                  //
//  Author: D.H. Wright (SLAC)                                            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
//

#include "G4HadronSim/LCLeptonPhysics.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// processes
//#include "G4MultipleScattering.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hIonisation.hh"
#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"

// models
#include "G4ElectroNuclearReaction.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoTau.hh"


LCLeptonPhysics::LCLeptonPhysics(const G4String& name, 
			   GlastMS::MultipleScatteringFactory& msFactory, GlastMS::EnergyLossFactory& eLossFactory): 
  G4VPhysicsConstructor(name), m_msFactory(msFactory), m_eLossFactory(eLossFactory)
{;}


LCLeptonPhysics::~LCLeptonPhysics()
{;}


void LCLeptonPhysics::ConstructParticle()
{ 
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4TauMinus::TauMinusDefinition();
  G4TauPlus::TauPlusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();
}


void LCLeptonPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;

  // Model for e+/e- nuclear reactions
   
  G4ElectroNuclearReaction* theElectronReaction =
                                   new G4ElectroNuclearReaction();

  // Electron physics

  pManager = G4Electron::Electron()->GetProcessManager();

  G4VContinuousDiscreteProcess* theElectronMultipleScattering = m_msFactory();
  G4VContinuousDiscreteProcess* theElectronIonisation         = 
    m_eLossFactory(GlastMS::EnergyLossFactory::ELECTRON, GlastMS::EnergyLossFactory::IONIZATION);
  G4VContinuousDiscreteProcess* theElectronBremsStrahlung     = 
    m_eLossFactory(GlastMS::EnergyLossFactory::ELECTRON, GlastMS::EnergyLossFactory::BREMSSTRAHLUNG);
  
  pManager->AddProcess(theElectronMultipleScattering, -1, 1, 1);
  pManager->AddProcess(theElectronIonisation,         -1, 2, 2);
  pManager->AddProcess(theElectronBremsStrahlung,     -1, 3, 3);

  //  pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
  //pManager->AddProcess(new G4eIonisation(),        -1, 2, 2);
  //pManager->AddProcess(new G4eBremsstrahlung(),    -1,-1, 3);  

  G4ElectronNuclearProcess* theElectronNuclearProcess =
                                   new G4ElectronNuclearProcess();
  theElectronNuclearProcess->RegisterMe(theElectronReaction);
  pManager->AddProcess(theElectronNuclearProcess, -1, -1, 4);

  // Positron physics

  pManager = G4Positron::Positron()->GetProcessManager(); 

  //  pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
  //pManager->AddProcess(new G4eIonisation(),        -1, 2, 2);
  //pManager->AddProcess(new G4eBremsstrahlung(),    -1,-1, 3);  
  //pManager->AddProcess(new G4eplusAnnihilation(),   0,-1, 4);
  
  G4VContinuousDiscreteProcess* thePositronMultipleScattering = m_msFactory();
  G4VContinuousDiscreteProcess* thePositronIonisation         = 
    m_eLossFactory(GlastMS::EnergyLossFactory::POSITRON, GlastMS::EnergyLossFactory::IONIZATION);
  G4VContinuousDiscreteProcess* thePositronBremsStrahlung     = 
    m_eLossFactory(GlastMS::EnergyLossFactory::POSITRON, GlastMS::EnergyLossFactory::BREMSSTRAHLUNG);
  G4eplusAnnihilation*          theAnnihilation               = new G4eplusAnnihilation();
  
  
  pManager->AddProcess(thePositronMultipleScattering, -1,  1, 1);
  pManager->AddProcess(thePositronIonisation,         -1,  2, 2);
  pManager->AddProcess(thePositronBremsStrahlung,     -1,  3, 3);
  pManager->AddProcess(theAnnihilation,                0, -1, 4);
  G4PositronNuclearProcess* thePositronNuclearProcess =
                                   new G4PositronNuclearProcess();
  thePositronNuclearProcess->RegisterMe(theElectronReaction);
  pManager->AddProcess(thePositronNuclearProcess, -1, -1, 5);

  // Muon-

  pManager = G4MuonMinus::MuonMinus()->GetProcessManager(); 

  G4VContinuousDiscreteProcess* fMuMinusIonisation         = 
    m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::IONIZATION);
  G4VContinuousDiscreteProcess* fMuMinusMultipleScattering = m_msFactory(); 
  G4VContinuousDiscreteProcess* fMuMinusBremsstrahlung     = 
    m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::BREMSSTRAHLUNG);
  G4VContinuousDiscreteProcess* fMuMinusPairProduction     = 
    m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::PAIRPRODUCTION);
  
  pManager->AddProcess(fMuMinusMultipleScattering, -1,           1, 1);
  pManager->AddProcess(fMuMinusIonisation,         -1,           2, 2);
  pManager->AddProcess(fMuMinusBremsstrahlung,     -1,           -1, 3);
  pManager->AddProcess(fMuMinusPairProduction,     -1,           -1, 4); 
  
  //  pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
  //pManager->AddProcess(new G4MuIonisation(),       -1, 2, 2);
  //pManager->AddProcess(new G4MuBremsstrahlung(),   -1,-1, 3);  
  //pManager->AddProcess(new G4MuPairProduction(),   -1,-1, 4);
  
  // Muon+

  pManager = G4MuonPlus::MuonPlus()->GetProcessManager(); 
  //  pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
  //pManager->AddProcess(new G4MuIonisation(),       -1, 2, 2);
  //pManager->AddProcess(new G4MuBremsstrahlung(),   -1,-1, 3);  
  //pManager->AddProcess(new G4MuPairProduction(),   -1,-1, 4);
  G4VContinuousDiscreteProcess* fMuPlusIonisation         = 
    m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::IONIZATION);
  G4VContinuousDiscreteProcess* fMuPlusMultipleScattering = m_msFactory(); 
  G4VContinuousDiscreteProcess* fMuPlusBremsstrahlung     = 
    m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::BREMSSTRAHLUNG);
  G4VContinuousDiscreteProcess* fMuPlusPairProduction     = 
    m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::PAIRPRODUCTION);
  
  pManager->AddProcess(fMuPlusMultipleScattering, -1,           1, 1);
  pManager->AddProcess(fMuPlusIonisation,         -1,           2, 2);
  pManager->AddProcess(fMuPlusBremsstrahlung,     -1,           -1, 3);
  pManager->AddProcess(fMuPlusPairProduction,     -1,           -1, 4); 
  

  // Tau-

  pManager = G4TauMinus::TauMinus()->GetProcessManager();
  G4VContinuousDiscreteProcess* theTauMinusMult = m_msFactory();
  //  pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
  pManager->AddProcess(theTauMinusMult, -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);
  
  // Tau+
  
  pManager = G4TauPlus::TauPlus()->GetProcessManager();

  G4VContinuousDiscreteProcess* theTauPlusMult = m_msFactory();
  //  pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
  pManager->AddProcess(theTauPlusMult, -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);


}
