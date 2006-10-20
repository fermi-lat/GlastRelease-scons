//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
// GEANT4 tag $Name$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmStandardPhysics72
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 19.12.2005 V.Ivanchenko rename 71 -> 72
// 15.06.2006 V.Ivanchenko use this class as a constructor of fast EM physics
//
//----------------------------------------------------------------------------
//

#include "G4EmStandardPhysics72.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering71.hh"
#include "G4MultipleScattering.hh"
#include "G4UrbanMscModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

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

G4EmStandardPhysics72::G4EmStandardPhysics72(const G4String& name, G4int ver,
   G4bool msc): G4VPhysicsConstructor(name), verbose(ver), mscStepLimit(msc)
{
  G4LossTableManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysics72::~G4EmStandardPhysics72()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysics72::ConstructParticle()
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

void G4EmStandardPhysics72::ConstructProcess()
{
  // Add standard EM Processes
  G4MultipleScattering* msc = 0;

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {

      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);

    } else if (particleName == "e-") {

      msc = new G4MultipleScattering();
      msc->MscStepLimitation(true,0.1);
      G4eIonisation* eion = new G4eIonisation;
      eion->ActivateSubCutoff(true);
      G4eBremsstrahlung* brem = new G4eBremsstrahlung();
      brem->ActivateSubCutoff(true);
      if(verbose > 1)
        G4cout << "### EmStandard72 instantiates eIoni and msc71 for " 
               << particleName << G4endl;
      pmanager->AddProcess(msc, -1, 1, 1);
      pmanager->AddProcess(eion,-1, 2, 2);
      pmanager->AddProcess(brem,-1,-3, 3);

    } else if (particleName == "e+") {

      msc = new G4MultipleScattering();
      msc->MscStepLimitation(true,0.1);
      if(verbose > 1)
        G4cout << "### EmStandard72 instantiates eIoni and msc71 for " 
               << particleName << G4endl;
      pmanager->AddProcess(msc, -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,          -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,      -1,-3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation,     0,-1, 4);

    } else if (particleName == "mu+" ||
               particleName == "mu-"    ) {

      msc = new G4MultipleScattering();
      msc->MscStepLimitation(false,0.2);
      if(verbose > 1)
        G4cout << "### EmStandard72 instantiates muIoni and msc71 for " 
               << particleName << G4endl;
      pmanager->AddProcess(msc,-1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation,        -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung,    -1,-3, 3);
      pmanager->AddProcess(new G4MuPairProduction,    -1,-4, 4);

    } else if (particleName == "alpha" ||
               particleName == "He3" ||
               particleName == "GenericIon") {

      msc = new G4MultipleScattering();
      msc->MscStepLimitation(false,0.2);
      if(verbose > 1)
        G4cout << "### EmStandard72 instantiates ionIoni and msc71 for " 
               << particleName << G4endl;
      pmanager->AddProcess(msc, -1, 1, 1);
      pmanager->AddProcess(new G4ionIonisation,        -1, 2, 2);

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
      if(verbose > 1)
        G4cout << "### EmStandard72 instantiates hIoni and msc71 for " 
               << particleName << G4endl;
      msc = new G4MultipleScattering();
      msc->MscStepLimitation(false,0.2);
      pmanager->AddProcess(msc, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,          -1, 2, 2);
    }
  }
  G4EmProcessOptions opt;
  opt.SetVerbose(verbose);
  //  opt.SetMscStepLimitation(false, 0.2);
  //  opt.SetSubCutoff(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
