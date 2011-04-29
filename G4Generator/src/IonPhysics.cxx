// File and Version Information:
// $Header$
//
// Description: This class manages the building of ions and their processes
//
// Author(s):
//      F.Longo & F.Paladin

#include "IonPhysics.h"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip> 

IonPhysics::IonPhysics(const G4String& name, std::string& physicsChoice
                       , GlastMS::MultipleScatteringFactory& msfactory)
  :  G4VPhysicsConstructor(name), m_physicsChoice(physicsChoice), m_msFactory(msfactory)
{
}

IonPhysics::~IonPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

// Nuclei
#include "G4IonConstructor.hh"

void IonPhysics::ConstructParticle()
{
  // Purpose and Method: this method is invoked by G4 to build the particles
  //                     classes

  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}


#include "G4ProcessManager.hh"
#include "G4HadronElasticProcess.hh"
#include "G4LElastic.hh"

#include "G4DeuteronInelasticProcess.hh"
#include "G4LEDeuteronInelastic.hh"

#include "G4TritonInelasticProcess.hh"
#include "G4LETritonInelastic.hh"

#include "G4AlphaInelasticProcess.hh"
#include "G4LEAlphaInelastic.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4hMultipleScattering.hh"

//--------------------------------------
// Ions hadronics models
//--------------------------------------

#include "G4IonInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "EpaxInelasticModel.hh"
#include "G4IonsShenCrossSection.hh"

void IonPhysics::ConstructProcess()
{
  // Purpose and Method: this method is invoked by G4 to build the physics
  //                     processes table

  G4ProcessManager * pManager = 0;
  
  // Full physics
  
  
  //  if (m_physicsChoice=="full")
  //  {
  // -----------------------
  // Ions Physics
  // ------------------------
  
  G4double theMin = 0;
  G4double theMax = 100*TeV; // corretto?? 
  G4IonsShenCrossSection* theIonXSec = new G4IonsShenCrossSection();
  
    // Elastic Process
  
  G4LElastic*            theIonElasticModel = new G4LElastic();
  
  // Generic Ion physics
  
  G4HadronElasticProcess* theIonElasticProcess = new G4HadronElasticProcess();
  theIonElasticProcess->RegisterMe(theIonElasticModel);
  G4IonInelasticProcess*  theIonInelasticProcess = new G4IonInelasticProcess();
  //   G4BinaryLightIonReaction* theIonModel = new  G4BinaryLightIonReaction();  
  EpaxInelasticModel* theIonModel = new  EpaxInelasticModel();
  theIonModel->SetMinEnergy(theMin);
  theIonModel->SetMaxEnergy(theMax);
  theIonInelasticProcess->RegisterMe(theIonModel);
  theIonInelasticProcess->AddDataSet(theIonXSec);  
  
    // Deuteron physics
  
  G4HadronElasticProcess*      theDElasticProcess = new G4HadronElasticProcess();
  theDElasticProcess->RegisterMe(theIonElasticModel);
  G4DeuteronInelasticProcess*  fDeuteronProcess = new  G4DeuteronInelasticProcess();
  G4LEDeuteronInelastic*      fDeuteronModel = new G4LEDeuteronInelastic();
  fDeuteronModel->SetMinEnergy(theMin);
  fDeuteronModel->SetMaxEnergy(theMax);
  fDeuteronProcess->RegisterMe(fDeuteronModel);
  fDeuteronProcess->AddDataSet(theIonXSec);  
  
    // Triton physics
  
  G4HadronElasticProcess*      theTElasticProcess = new  G4HadronElasticProcess();
  theTElasticProcess->RegisterMe(theIonElasticModel);
  G4TritonInelasticProcess*    fTritonProcess = new G4TritonInelasticProcess();
  G4LETritonInelastic*        fTritonModel= new G4LETritonInelastic();
  fTritonModel->SetMinEnergy(theMin);
  fTritonModel->SetMaxEnergy(theMax);
  fTritonProcess->RegisterMe(fTritonModel);
  fTritonProcess->AddDataSet(theIonXSec);    
  
  // Alpha physics
  
  G4HadronElasticProcess*      theAElasticProcess = new G4HadronElasticProcess();
  theAElasticProcess->RegisterMe(theIonElasticModel);
  G4AlphaInelasticProcess*     fAlphaProcess = new G4AlphaInelasticProcess();
  G4LEAlphaInelastic*         fAlphaModel = new G4LEAlphaInelastic();
  fAlphaModel->SetMinEnergy(theMin);
  fAlphaModel->SetMaxEnergy(theMax);
  fAlphaProcess->RegisterMe(fAlphaModel);
  fAlphaProcess->AddDataSet(theIonXSec);  
  
  // He3 physics
    
  G4HadronElasticProcess* theHe3ElasticProcess = new G4HadronElasticProcess();
  theHe3ElasticProcess->RegisterMe(theIonElasticModel);
  
  
  // Building processes
   
  // Generic Ion
  
  pManager = G4GenericIon::GenericIon()->GetProcessManager();
  //pManager->AddDiscreteProcess(theIonElasticProcess);
  pManager->AddDiscreteProcess(theIonInelasticProcess);
  
  // Deuteron 
  
  pManager = G4Deuteron::Deuteron()->GetProcessManager();
  pManager->AddDiscreteProcess(theDElasticProcess);
  pManager->AddDiscreteProcess(fDeuteronProcess);
  
  // Triton 
  pManager = G4Triton::Triton()->GetProcessManager();
  pManager->AddDiscreteProcess(theTElasticProcess);
  pManager->AddDiscreteProcess(fTritonProcess);
  
  // Alpha 
  pManager = G4Alpha::Alpha()->GetProcessManager();
  pManager->AddDiscreteProcess(theAElasticProcess);
  pManager->AddDiscreteProcess(fAlphaProcess);
  
  // He3
  pManager = G4He3::He3()->GetProcessManager();
  pManager->AddDiscreteProcess(theHe3ElasticProcess);


  //    }

  
  // Generic Ion 
  
  pManager = G4GenericIon::GenericIon()->GetProcessManager();
  
  // G4VContinuousDiscreteProcess*   fIonMultipleScattering = m_msFactory();
  G4ionIonisation*          fIonIonisation = new  G4ionIonisation(); 
  //G4hIonisation*          fIonIonisation = new  G4hIonisation();
  //G4hIonisation52*          fIonIonisation = new  G4hIonisation52();
  // G4hLowEnergyIonisation* fIonIonisation = new G4hLowEnergyIonisation();
  pManager->AddProcess(fIonIonisation, ordInActive, 2, 2);
  G4VProcess* fIonMultipleScattering = new G4hMultipleScattering();
  pManager->AddProcess(fIonMultipleScattering);
  pManager->SetProcessOrdering(fIonMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fIonMultipleScattering, idxPostStep,  1);

  // Deuteron 
  
  pManager = G4Deuteron::Deuteron()->GetProcessManager();

  G4VContinuousDiscreteProcess*        fDeuteronMultipleScattering = m_msFactory();
  G4hIonisation*               fDeuteronIonisation= new  G4hIonisation();
  pManager->AddProcess(fDeuteronIonisation, ordInActive, 2, 2);
  pManager->AddProcess(fDeuteronMultipleScattering);
  pManager->SetProcessOrdering(fDeuteronMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fDeuteronMultipleScattering, idxPostStep,  1);

  // Triton
      
  pManager = G4Triton::Triton()->GetProcessManager();

  G4VContinuousDiscreteProcess*        fTritonMultipleScattering = m_msFactory();
  G4hIonisation*               fTritonIonisation = new  G4hIonisation();
  pManager->AddProcess(fTritonIonisation, ordInActive, 2, 2);
  pManager->AddProcess(fTritonMultipleScattering);
  pManager->SetProcessOrdering(fTritonMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fTritonMultipleScattering, idxPostStep,  1);  

  // Alpha 
  
  pManager = G4Alpha::Alpha()->GetProcessManager();

  G4VContinuousDiscreteProcess*        fAlphaMultipleScattering  = m_msFactory();
  G4hIonisation*               fAlphaIonisation= new  G4hIonisation();
  pManager->AddProcess(fAlphaIonisation, ordInActive, 2, 2);
  pManager->AddProcess(fAlphaMultipleScattering);
  pManager->SetProcessOrdering(fAlphaMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fAlphaMultipleScattering, idxPostStep,  1);

  // He3
  
  pManager = G4He3::He3()->GetProcessManager();
  G4VContinuousDiscreteProcess*        fHe3MultipleScattering  = m_msFactory();
  G4hIonisation*               fHe3Ionisation= new  G4hIonisation();
  pManager->AddProcess(fHe3Ionisation, ordInActive, 2, 2);
  pManager->AddProcess(fHe3MultipleScattering);
  pManager->SetProcessOrdering(fHe3MultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fHe3MultipleScattering, idxPostStep,  1); 

}      









