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


EMPhysics::EMPhysics(const G4String& name, 
                     GlastMS::MultipleScatteringFactory& msFactory,
                     GlastMS::EnergyLossFactory& eLossFactory )
               :  G4VPhysicsConstructor(name), m_msFactory(msFactory), m_eLossFactory(eLossFactory)
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
#include "G4eplusAnnihilation.hh"

#include "G4eBremsstrahlungModel.hh"
#include "G4VEnergyLossProcess.hh"

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
  G4VContinuousDiscreteProcess* theElectronIonisation         = 
      m_eLossFactory(GlastMS::EnergyLossFactory::ELECTRON, GlastMS::EnergyLossFactory::IONIZATION);
  G4VContinuousDiscreteProcess* theElectronBremsStrahlung     = 
      m_eLossFactory(GlastMS::EnergyLossFactory::ELECTRON, GlastMS::EnergyLossFactory::BREMSSTRAHLUNG);

    G4eBremsstrahlungModel* bm = new G4eBremsstrahlungModel();
    //bm->SetLPMflag(false);

    G4VEnergyLossProcess* ebrem = 
      reinterpret_cast<G4VEnergyLossProcess*>(theElectronBremsStrahlung);
    ebrem->AddEmModel(0,bm);

  pManager->AddProcess(theElectronMultipleScattering, -1, 1, 1);
  pManager->AddProcess(theElectronIonisation,         -1, 2, 2);
  pManager->AddProcess(theElectronBremsStrahlung,     -1, 3, 3);

  //Positron Physics

  pManager = G4Positron::Positron()->GetProcessManager();

  G4VContinuousDiscreteProcess* thePositronMultipleScattering = m_msFactory();
  G4VContinuousDiscreteProcess* thePositronIonisation         = 
      m_eLossFactory(GlastMS::EnergyLossFactory::POSITRON, GlastMS::EnergyLossFactory::IONIZATION);
  G4VContinuousDiscreteProcess* thePositronBremsStrahlung     = 
      m_eLossFactory(GlastMS::EnergyLossFactory::POSITRON, GlastMS::EnergyLossFactory::BREMSSTRAHLUNG);
  G4eplusAnnihilation*          theAnnihilation               = new G4eplusAnnihilation();

    G4eBremsstrahlungModel* bp = new G4eBremsstrahlungModel();
    //bp->SetLPMflag(false);
    G4VEnergyLossProcess* epbrem =
      reinterpret_cast<G4VEnergyLossProcess*>(thePositronBremsStrahlung);
    epbrem->AddEmModel(0,bp);

  pManager->AddProcess(thePositronMultipleScattering, -1,  1, 1);
  pManager->AddProcess(thePositronIonisation,         -1,  2, 2);
  pManager->AddProcess(thePositronBremsStrahlung,     -1,  3, 3);
  pManager->AddProcess(theAnnihilation,                0, -1, 4);
  

}



