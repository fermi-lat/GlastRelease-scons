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

#include "G4ionIonisation.hh"
#include "G4hIonisation.hh"
//THB#include "G4MultipleScattering.hh"


void IonPhysics::ConstructProcess()
{
  // Purpose and Method: this method is invoked by G4 to build the physics
  //                     processes table

  G4ProcessManager * pManager = 0;
  
  // Full physics
  
  if (m_physicsChoice=="full")
    {
      G4HadronElasticProcess* theElasticProcess = new  G4HadronElasticProcess();
      G4LElastic*            theElasticModel;
      theElasticModel = new G4LElastic();
      theElasticProcess->RegisterMe(theElasticModel);    

      // Generic Ion 
      
      pManager = G4GenericIon::GenericIon()->GetProcessManager();

      pManager->AddDiscreteProcess(theElasticProcess);

      // Deuteron 

      pManager = G4Deuteron::Deuteron()->GetProcessManager();

      pManager->AddDiscreteProcess(theElasticProcess);
      G4DeuteronInelasticProcess*  fDeuteronProcess = new  G4DeuteronInelasticProcess();
      G4LEDeuteronInelastic*      fDeuteronModel;
      fDeuteronModel = new G4LEDeuteronInelastic();
      fDeuteronProcess->RegisterMe(fDeuteronModel);
      pManager->AddDiscreteProcess(fDeuteronProcess);

      // Triton
      
      pManager = G4Triton::Triton()->GetProcessManager();

      pManager->AddDiscreteProcess(theElasticProcess);
      G4TritonInelasticProcess*    fTritonProcess = new G4TritonInelasticProcess();
      G4LETritonInelastic*        fTritonModel;
      fTritonModel = new G4LETritonInelastic();
      fTritonProcess->RegisterMe(fTritonModel);
      pManager->AddDiscreteProcess(fTritonProcess);

      // Alpha 
      
      pManager = G4Alpha::Alpha()->GetProcessManager();
 
      pManager->AddDiscreteProcess(theElasticProcess);
      G4AlphaInelasticProcess*     fAlphaProcess = new G4AlphaInelasticProcess();
      G4LEAlphaInelastic*         fAlphaModel;   
      fAlphaModel = new G4LEAlphaInelastic();
      fAlphaProcess->RegisterMe(fAlphaModel);
      pManager->AddDiscreteProcess(fAlphaProcess);

      // He3
      
      pManager = G4He3::He3()->GetProcessManager();
      
      pManager->AddDiscreteProcess(theElasticProcess);

    }

  
  // Generic Ion 
  
  pManager = G4GenericIon::GenericIon()->GetProcessManager();
  
  G4VContinuousDiscreteProcess*   fIonMultipleScattering = m_msFactory();
  G4ionIonisation*          fIonIonisation = new  G4ionIonisation();
  pManager->AddProcess(fIonIonisation, ordInActive, 2, 2);
  
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
  G4ionIonisation*               fAlphaIonisation= new  G4ionIonisation();
  pManager->AddProcess(fAlphaIonisation, ordInActive, 2, 2);
  pManager->AddProcess(fAlphaMultipleScattering);
  pManager->SetProcessOrdering(fAlphaMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fAlphaMultipleScattering, idxPostStep,  1);

  // He3
  
  pManager = G4He3::He3()->GetProcessManager();
  G4VContinuousDiscreteProcess*        fHe3MultipleScattering  = m_msFactory();
  G4ionIonisation*               fHe3Ionisation= new  G4ionIonisation();
  pManager->AddProcess(fHe3Ionisation, ordInActive, 2, 2);
  pManager->AddProcess(fHe3MultipleScattering);
  pManager->SetProcessOrdering(fHe3MultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fHe3MultipleScattering, idxPostStep,  1); 

}










