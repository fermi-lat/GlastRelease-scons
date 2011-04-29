//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id$
// GEANT4 tag $Name$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmPenelopePhysics
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 05.12.2005 V.Ivanchenko add controlled verbosity
//
//----------------------------------------------------------------------------
//

#include "G4HadronSim/G4EmPenelopePhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh"

#include "G4PenelopeCompton.hh"
#include "G4PenelopeGammaConversion.hh"
#include "G4PenelopePhotoElectric.hh"
#include "G4PenelopeRayleigh.hh"

//#include "G4MultipleScattering.hh"

#include "G4eMultipleScattering.hh"
#include "G4PenelopeIonisation.hh"
#include "G4PenelopeBremsstrahlung.hh"
#include "G4PenelopeAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmPenelopePhysics::G4EmPenelopePhysics(const G4String& name, G4int ver,
					 //   G4bool msc,
					 GlastMS::MultipleScatteringFactory& msFactory,
					 GlastMS::EnergyLossFactory& eLossFactory ): 
  G4VPhysicsConstructor(name), verbose(ver), m_msFactory(msFactory), m_eLossFactory(eLossFactory)
  //mscStepLimit(msc),
{
  G4LossTableManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmPenelopePhysics::~G4EmPenelopePhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmPenelopePhysics::ConstructParticle()
{
// gamma
  G4Gamma::Gamma();

// leptons
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();

// mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();

// barions
  G4Proton::Proton();
  G4AntiProton::AntiProton();

// ions
  G4Deuteron::Deuteron();
  G4Triton::Triton();
  G4He3::He3();
  G4Alpha::Alpha();
  G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmPenelopePhysics::ConstructProcess()
{
  // Add standard EM Processes

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {

      pmanager -> AddDiscreteProcess(new G4PenelopePhotoElectric);
      pmanager -> AddDiscreteProcess(new G4PenelopeCompton);
      pmanager -> AddDiscreteProcess(new G4PenelopeGammaConversion);
      pmanager -> AddDiscreteProcess(new G4PenelopeRayleigh);

    } else if (particleName == "e-") {

      //      if(verbose > 1)
      // G4cout << "### EmPenelope instantiates eIoni and msc80 for " 
      //         << particleName << G4endl;
      //pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
      //pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
      //pmanager->AddProcess(new G4eBremsstrahlung(),  -1, 3, 3);
      G4VContinuousDiscreteProcess* theElectronMultipleScattering = m_msFactory();
      // G4VContinuousDiscreteProcess* theElectronIonisation         = 
      //	m_eLossFactory(GlastMS::EnergyLossFactory::ELECTRON, GlastMS::EnergyLossFactory::IONIZATION);
      //G4VContinuousDiscreteProcess* theElectronBremsStrahlung     = 
      //m_eLossFactory(GlastMS::EnergyLossFactory::ELECTRON, GlastMS::EnergyLossFactory::BREMSSTRAHLUNG);
      
      pmanager->AddProcess(theElectronMultipleScattering, -1, 1, 1);

      G4PenelopeIonisation* ioni = new G4PenelopeIonisation();
      G4PenelopeBremsstrahlung* brem = new G4PenelopeBremsstrahlung();
      pmanager -> AddProcess(ioni, -1, 2,2);
      pmanager -> AddProcess(brem, -1,-1,3);      
      
    } else if (particleName == "e+") {

      /*    if(verbose > 1)
        G4cout << "### EmPenelope instantiates eIoni and msc80 for " 
	<< particleName << G4endl;
	pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
	pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
	pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3, 3);
	pmanager->AddProcess(new G4eplusAnnihilation,   0,-1, 4);
      */
      G4VContinuousDiscreteProcess* thePositronMultipleScattering = m_msFactory();
      //G4VContinuousDiscreteProcess* thePositronIonisation         = 
      //m_eLossFactory(GlastMS::EnergyLossFactory::POSITRON, GlastMS::EnergyLossFactory::IONIZATION);
      //G4VContinuousDiscreteProcess* thePositronBremsStrahlung     = 
      //m_eLossFactory(GlastMS::EnergyLossFactory::POSITRON, GlastMS::EnergyLossFactory::BREMSSTRAHLUNG);
      //G4eplusAnnihilation*          theAnnihilation               = new G4eplusAnnihilation();
      
      pmanager->AddProcess(thePositronMultipleScattering, -1,  1, 1);
      //pmanager->AddProcess(thePositronIonisation,         -1,  2, 2);
      //pmanager->AddProcess(thePositronBremsStrahlung,     -1,  3, 3);
      //pmanager->AddProcess(theAnnihilation,                0, -1, 4);
      pmanager -> AddProcess(new G4PenelopeIonisation,    -1, 2,2);
      pmanager -> AddProcess(new G4PenelopeBremsstrahlung,-1,-1,3);
      pmanager -> AddProcess(new G4PenelopeAnnihilation,   0,-1,4);
      

    } else if (particleName == "mu+")
      //               particleName == "mu-") 
    {

      /*     if(verbose > 1)
        G4cout << "### EmPenelope instantiates muIoni and msc80 for " 
	<< particleName << G4endl;
	pmanager->AddProcess(new G4MultipleScattering,-1, 1, 1);
	pmanager->AddProcess(new G4MuIonisation,      -1, 2, 2);
	pmanager->AddProcess(new G4MuBremsstrahlung,  -1, 3, 3);
	pmanager->AddProcess(new G4MuPairProduction,  -1, 4, 4);
      */
      G4VContinuousDiscreteProcess* fMuPlusIonisation         = 
	m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::IONIZATION);
      G4VContinuousDiscreteProcess* fMuPlusMultipleScattering = m_msFactory(); 
      G4VContinuousDiscreteProcess* fMuPlusBremsstrahlung     = 
	m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::BREMSSTRAHLUNG);
      G4VContinuousDiscreteProcess* fMuPlusPairProduction     = 
	m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::PAIRPRODUCTION);
            
      pmanager->AddProcess(fMuPlusMultipleScattering, -1,           1, 1);
      pmanager->AddProcess(fMuPlusIonisation,         -1,           2, 2);
      pmanager->AddProcess(fMuPlusBremsstrahlung,     -1,           3, 3);
      pmanager->AddProcess(fMuPlusPairProduction,     -1,           4, 4);      
    }
    else if (particleName == "mu-")
      {
	
	/*     if(verbose > 1)
	       G4cout << "### EmPenelope instantiates muIoni and msc80 for " 
	       << particleName << G4endl;
	       pmanager->AddProcess(new G4MultipleScattering,-1, 1, 1);
	       pmanager->AddProcess(new G4MuIonisation,      -1, 2, 2);
	       pmanager->AddProcess(new G4MuBremsstrahlung,  -1, 3, 3);
	       pmanager->AddProcess(new G4MuPairProduction,  -1, 4, 4);
	*/
	G4VContinuousDiscreteProcess* fMuMinusIonisation         = 
	  m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::IONIZATION);
	G4VContinuousDiscreteProcess* fMuMinusMultipleScattering = m_msFactory(); 
	G4VContinuousDiscreteProcess* fMuMinusBremsstrahlung     = 
	  m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::BREMSSTRAHLUNG);
	G4VContinuousDiscreteProcess* fMuMinusPairProduction     = 
	  m_eLossFactory(GlastMS::EnergyLossFactory::MUON, GlastMS::EnergyLossFactory::PAIRPRODUCTION);
		
	pmanager->AddProcess(fMuMinusMultipleScattering, -1,           1, 1);
	pmanager->AddProcess(fMuMinusIonisation,         -1,           2, 2);
	pmanager->AddProcess(fMuMinusBremsstrahlung,     -1,           3, 3);     
	pmanager->AddProcess(fMuMinusPairProduction,     -1,           4, 4);

      }
    else if (particleName == "alpha" ||
	     particleName == "He3" ||
	     particleName == "GenericIon") {
      
      //      if(verbose > 1)
      //  G4cout << "### EmPenelope instantiates ionIoni and msc80 for " 
      //         << particleName << G4endl;
      //      pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
      G4VContinuousDiscreteProcess* theIonMult = m_msFactory();
      pmanager->AddProcess(theIonMult, -1, 1, 1);
      pmanager->AddProcess(new G4ionIonisation,      -1, 2, 2);
      
    } else if (particleName == "anti_omega-" ||
               particleName == "anti_proton" ||
               particleName == "anti_sigma+" ||
               particleName == "anti_sigma-" ||
               particleName == "anti_xi-" ||
               particleName == "deuteron" ||
               particleName == "kaon+" ||
               particleName == "kaon-" ||
               particleName == "omega-" ||
               particleName == "pi+" ||
               particleName == "pi-" ||
               particleName == "proton" ||
               particleName == "sigma+" ||
               particleName == "sigma-" ||
               particleName == "tau+" ||
               particleName == "tau-" ||
               particleName == "triton" ||
               particleName == "xi-" ) {

      //if(verbose > 1)
      //  G4cout << "### EmPenelope instantiates hIoni and msc80 for " 
      //         << particleName << G4endl;
      //pmanager->AddProcess(new G4MultipleScattering,-1, 1, 1);
      G4VContinuousDiscreteProcess* theHadronMult = m_msFactory();
      pmanager->AddProcess(theHadronMult,-1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,       -1, 2, 2);
    }
  }
  /*
  G4EmProcessOptions opt; // to be checked
  opt.SetVerbose(verbose);
  if(!mscStepLimit) 
  opt.SetMscStepLimitation(false);
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
