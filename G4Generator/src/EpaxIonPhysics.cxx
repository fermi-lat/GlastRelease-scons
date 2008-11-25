// File and Version Information:
// $Header$
//
// Description: This class manages the building of ions and their processes
//
// Author(s):
//      F.Longo & F.Paladin

#include "EpaxIonPhysics.h"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip> 

EpaxIonPhysics::EpaxIonPhysics(const G4String& name, std::string& physicsChoice
                       , GlastMS::MultipleScatteringFactory& msfactory)
  :  G4VPhysicsConstructor(name), m_physicsChoice(physicsChoice), m_msFactory(msfactory)
{
}

EpaxIonPhysics::~EpaxIonPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

// Nuclei
#include "G4IonConstructor.hh"

void EpaxIonPhysics::ConstructParticle()
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

//--------------------------------------
// Ions hadronics models
//--------------------------------------

#include "G4IonInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "EpaxInelasticModel.hh"
#include "G4IonsShenCrossSection.hh"

void EpaxIonPhysics::ConstructProcess()
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
  //G4BinaryLightIonReaction* theIonModel = new  G4BinaryLightIonReaction();  
  EpaxInelasticModel* theIonModel = new  EpaxInelasticModel();
  //theIonModel->SetMinEnergy(theMin);
  //theIonModel->SetMaxEnergy(theMax);
  theIonInelasticProcess->RegisterMe(theIonModel);
  theIonInelasticProcess->AddDataSet(theIonXSec);  
  
    // Deuteron physics
  
  G4HadronElasticProcess*      theDElasticProcess = new G4HadronElasticProcess();
  theDElasticProcess->RegisterMe(theIonElasticModel);
  G4DeuteronInelasticProcess*  fDeuteronProcess = new  G4DeuteronInelasticProcess();
  G4LEDeuteronInelastic*      fDeuteronModel = new G4LEDeuteronInelastic();
  //  fDeuteronModel->SetMinEnergy(theMin);
  //fDeuteronModel->SetMaxEnergy(theMax);
  fDeuteronProcess->RegisterMe(fDeuteronModel);
  fDeuteronProcess->AddDataSet(theIonXSec);  
  
    // Triton physics
  
  G4HadronElasticProcess*      theTElasticProcess = new  G4HadronElasticProcess();
  theTElasticProcess->RegisterMe(theIonElasticModel);
  G4TritonInelasticProcess*    fTritonProcess = new G4TritonInelasticProcess();
  G4LETritonInelastic*        fTritonModel= new G4LETritonInelastic();
  //fTritonModel->SetMinEnergy(theMin);
  //fTritonModel->SetMaxEnergy(theMax);
  fTritonProcess->RegisterMe(fTritonModel);
  fTritonProcess->AddDataSet(theIonXSec);    
  
  // Alpha physics
  
  G4HadronElasticProcess*      theAElasticProcess = new G4HadronElasticProcess();
  theAElasticProcess->RegisterMe(theIonElasticModel);
  G4AlphaInelasticProcess*     fAlphaProcess = new G4AlphaInelasticProcess();
  G4LEAlphaInelastic*         fAlphaModel = new G4LEAlphaInelastic();
  //EpaxInelasticModel* theAlphaModel = new  EpaxInelasticModel();
  //theIonModel->SetMinEnergy(theMin);
  //theIonModel->SetMaxEnergy(theMax);
  //fAlphaProcess->RegisterMe(theAlphaModel);
  fAlphaProcess->RegisterMe(fAlphaModel);
  fAlphaProcess->AddDataSet(theIonXSec);  
  
  // He3 physics
    
  G4HadronElasticProcess* theHe3ElasticProcess = new G4HadronElasticProcess();
  theHe3ElasticProcess->RegisterMe(theIonElasticModel);
  
  
  // Building processes
   
  // Generic Ion
  
  pManager = G4GenericIon::GenericIon()->GetProcessManager();
  pManager->AddDiscreteProcess(theIonElasticProcess);
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
      
}      









