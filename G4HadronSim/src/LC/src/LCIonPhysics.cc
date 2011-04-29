//
////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Title:  Ion physics for a Linear Collider Detector                    //
//  Date:   7 July 2004                                                   //
//  Author: D.H. Wright (SLAC)                                            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
//

#include "G4HadronSim/LCIonPhysics.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4IonConstructor.hh"

// processes

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4HadronElasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

// models
#include "G4LElastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"

LCIonPhysics::LCIonPhysics(const G4String& name, 
			   GlastMS::MultipleScatteringFactory& msFactory)
  :  G4VPhysicsConstructor(name), m_msFactory(msFactory)
{;}


LCIonPhysics::~LCIonPhysics()
{;}


void LCIonPhysics::ConstructParticle()
{ 
  // Construct light ions (d, t, 3He, alpha, and generic ion)
  G4IonConstructor ionConstruct;
  ionConstruct.ConstructParticle();
}


void LCIonPhysics::ConstructProcess()
{
  // Hadronic Elastic Process and Model (for all ions except generic ion)

  G4HadronElasticProcess* elasticProcess = new G4HadronElasticProcess();
  G4LElastic* elasticModel = new G4LElastic();
  elasticProcess->RegisterMe(elasticModel);

  // Hadronic inelastic models

  G4ProcessManager * pManager = 0;
  
  ///////////////////
  //               //
  //   Deuteron    //
  //               //
  ///////////////////

  pManager = G4Deuteron::Deuteron()->GetProcessManager();

  // EM processes

  G4VContinuousDiscreteProcess*   theDeuteronMult = m_msFactory();
  //  pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
  pManager->AddProcess(theDeuteronMult, -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4DeuteronInelasticProcess* dinelProc = new G4DeuteronInelasticProcess();
  G4LEDeuteronInelastic* LEPdModel = new G4LEDeuteronInelastic();
  dinelProc->RegisterMe(LEPdModel);
  pManager->AddDiscreteProcess(dinelProc);

  ///////////////////
  //               //
  //    Triton     //
  //               //
  ///////////////////

  pManager = G4Triton::Triton()->GetProcessManager(); 

  // EM processes
  G4VContinuousDiscreteProcess*   theTritonMult = m_msFactory();
  //  pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
  pManager->AddProcess(theTritonMult, -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);


  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4TritonInelasticProcess* tinelProc = new G4TritonInelasticProcess();
  G4LETritonInelastic* LEPtModel = new G4LETritonInelastic();
  tinelProc->RegisterMe(LEPtModel);
  pManager->AddDiscreteProcess(tinelProc);

  ///////////////////
  //               //
  //      3He      //
  //               //
  ///////////////////

  pManager = G4He3::He3()->GetProcessManager(); 

  // EM processes

  G4VContinuousDiscreteProcess*   theHe3Mult = m_msFactory();
  //  pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
  pManager->AddProcess(theHe3Mult, -1, 1, 1); // bug? 
  pManager->AddProcess(new G4ionIonisation(),        -1, 2, 2);
  
  // hadron elastic 
  pManager->AddDiscreteProcess(elasticProcess);

  // NO INELASTIC PROCESS AVAILABLE FOR 3HE

  ///////////////////
  //               //
  //     Alpha     //
  //               //
  ///////////////////

  pManager = G4Alpha::Alpha()->GetProcessManager(); 

  // EM processes

  G4VContinuousDiscreteProcess*   theAlphaMult = m_msFactory();
  //  pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
  pManager->AddProcess(theAlphaMult, -1, 1, 1); // bug? 
  pManager->AddProcess(new G4ionIonisation(),        -1, 2, 2);

   // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AlphaInelasticProcess* ainelProc = new G4AlphaInelasticProcess();
  G4LEAlphaInelastic* LEPaModel = new G4LEAlphaInelastic();
  ainelProc->RegisterMe(LEPaModel);
  pManager->AddDiscreteProcess(ainelProc);

  ///////////////////
  //               //
  //  generic ion  //
  //               //
  ///////////////////

  pManager = G4GenericIon::GenericIon()->GetProcessManager();

  // Only EM processes for generic ion

  G4VContinuousDiscreteProcess*   theIonMult = m_msFactory();
  //  pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
  pManager->AddProcess(theIonMult, -1, 1, 1); // bug? 
  pManager->AddProcess(new G4ionIonisation(),        -1, 2, 2);
 
}
