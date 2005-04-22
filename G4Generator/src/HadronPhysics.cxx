


// File and Version Information:
// $Header$
//
// Description: This class manages the building of hadrons and their processes
//
// Author(s):
//      F.Longo & F.Paladin

#include "HadronPhysics.h"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>

HadronPhysics::HadronPhysics(const G4String& name, std::string& physicsChoice,
                             GlastMS::MultipleScatteringFactory& msFactory)
  :  G4VPhysicsConstructor(name),  m_physicsChoice(physicsChoice), m_msFactory(msFactory)
{
}

HadronPhysics::~HadronPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

// Nuclei
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

void HadronPhysics::ConstructParticle()
{
  // Purpose and Method: this method is invoked by G4 to build the particles
  //                     classes

  //  Construct all mesons
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  //  Construct all barions
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  

}


#include "G4ProcessManager.hh"
//THB #include "G4MultipleScattering.hh"
#include "G4hIonisation.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4LambdaInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4XiZeroInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"

// Low-energy Models
#include "G4LElastic.hh"   
#include "G4LFission.hh"
#include "G4LCapture.hh"

#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4LELambdaInelastic.hh"
#include "G4LEAntiLambdaInelastic.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "G4LEAntiSigmaMinusInelastic.hh"
#include "G4LEXiZeroInelastic.hh"
#include "G4LEXiMinusInelastic.hh"
#include "G4LEAntiXiZeroInelastic.hh"
#include "G4LEAntiXiMinusInelastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4LEOmegaMinusInelastic.hh"
#include "G4LEAntiOmegaMinusInelastic.hh"

// High-energy Models

#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4HENeutronInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"
#include "G4HELambdaInelastic.hh"
#include "G4HEAntiLambdaInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"
#include "G4HEAntiSigmaPlusInelastic.hh"
#include "G4HEAntiSigmaMinusInelastic.hh"
#include "G4HEXiZeroInelastic.hh"
#include "G4HEXiMinusInelastic.hh"
#include "G4HEAntiXiZeroInelastic.hh"
#include "G4HEAntiXiMinusInelastic.hh"
#include "G4HEOmegaMinusInelastic.hh"
#include "G4HEAntiOmegaMinusInelastic.hh"

// Stopping processes
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"

#ifdef TRIUMF_STOP_PIMINUS
#include "G4PionMinusAbsorptionAtRest.hh"
#else
#include "G4PiMinusAbsorptionAtRest.hh"
#endif
#ifdef TRIUMF_STOP_KMINUS
#include "G4KaonMinusAbsorption.hh"
#else
#include "G4KaonMinusAbsorptionAtRest.hh"
#endif

void HadronPhysics::ConstructProcess()
{
  // Purpose and Method: this method is invoked by G4 to build the physics
  //                     processes table

  G4ProcessManager * pManager = 0;
  
  
  // Only EM

  // PionPlus
  
  pManager = G4PionPlus::PionPlus()->GetProcessManager();
  
  G4VContinuousDiscreteProcess* thePionPlusMult = m_msFactory();
  G4hIonisation* thePionPlusIonisation = new G4hIonisation();
  pManager->AddProcess(thePionPlusIonisation, ordInActive,2, 2);
  pManager->AddProcess(thePionPlusMult);
  pManager->SetProcessOrdering(thePionPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(thePionPlusMult, idxPostStep, 1);  

 // PionMinus

  pManager = G4PionMinus::PionMinus()->GetProcessManager();

  G4VContinuousDiscreteProcess* thePionMinusMult = m_msFactory();
  G4hIonisation* thePionMinusIonisation =  new G4hIonisation();
  pManager->AddProcess(thePionMinusIonisation, ordInActive,2, 2);
  pManager->AddProcess(thePionMinusMult);
  pManager->SetProcessOrdering(thePionMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(thePionMinusMult, idxPostStep, 1);

  // KaonPlus

  pManager = G4KaonPlus::KaonPlus()->GetProcessManager();

  G4VContinuousDiscreteProcess* theKaonPlusMult = m_msFactory();
  G4hIonisation* theKaonPlusIonisation = new G4hIonisation();
  pManager->AddProcess(theKaonPlusIonisation, ordInActive,2, 2);
  pManager->AddProcess(theKaonPlusMult);
  pManager->SetProcessOrdering(theKaonPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theKaonPlusMult, idxPostStep, 1);

  // KaonMinus
  
  pManager = G4KaonMinus::KaonMinus()->GetProcessManager();
  
  G4VContinuousDiscreteProcess* theKaonMinusMult = m_msFactory();
  G4hIonisation* theKaonMinusIonisation = new G4hIonisation(); 
  pManager->AddProcess(theKaonMinusIonisation, ordInActive,2, 2);
  pManager->AddProcess(theKaonMinusMult);
  pManager->SetProcessOrdering(theKaonMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theKaonMinusMult, idxPostStep, 1);

  // Proton

  pManager = G4Proton::Proton()->GetProcessManager();

  G4VContinuousDiscreteProcess* theProtonMult = m_msFactory();
  G4hIonisation* theProtonIonisation = new G4hIonisation() ;
  pManager->AddProcess(theProtonIonisation, ordInActive,2, 2);
  pManager->AddProcess(theProtonMult);
  pManager->SetProcessOrdering(theProtonMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theProtonMult, idxPostStep, 1);

  // AntiProton
  
  pManager = G4AntiProton::AntiProton()->GetProcessManager();
  
  G4VContinuousDiscreteProcess* theAntiProtonMult = m_msFactory();
  G4hIonisation* theAntiProtonIonisation = new G4hIonisation();
  G4AntiProtonAnnihilationAtRest*  theAntiProtonAnnihilation = new G4AntiProtonAnnihilationAtRest();
  pManager->AddProcess(theAntiProtonIonisation, ordInActive,2, 2);
  pManager->AddProcess(theAntiProtonMult);
  pManager->SetProcessOrdering(theAntiProtonMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theAntiProtonMult, idxPostStep, 1);
  pManager->AddRestProcess(theAntiProtonAnnihilation);

  // AntiNeutron
  
  pManager = G4AntiNeutron::AntiNeutron()->GetProcessManager();

  G4AntiNeutronAnnihilationAtRest*  theAntiNeutronAnnihilation = new G4AntiNeutronAnnihilationAtRest();
  pManager->AddRestProcess(theAntiNeutronAnnihilation);  
  
  // SigmaMinus
  
  pManager = G4SigmaMinus::SigmaMinus()->GetProcessManager();
  
  G4VContinuousDiscreteProcess* theSigmaMinusMult =  m_msFactory();
  G4hIonisation* theSigmaMinusIonisation = new G4hIonisation();
  pManager->AddProcess(theSigmaMinusIonisation, ordInActive,2, 2);
  pManager->AddProcess(theSigmaMinusMult);
  pManager->SetProcessOrdering(theSigmaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theSigmaMinusMult, idxPostStep, 1);
 

// AntiSigmaMinus
  
  pManager = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();

  G4VContinuousDiscreteProcess* theAntiSigmaMinusMult = m_msFactory();
  G4hIonisation* theAntiSigmaMinusIonisation = new G4hIonisation();
  pManager->AddProcess(theAntiSigmaMinusIonisation, ordInActive,2, 2);
  pManager->AddProcess(theAntiSigmaMinusMult);
  pManager->SetProcessOrdering(theAntiSigmaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theAntiSigmaMinusMult, idxPostStep, 1);

  // SigmaPlus
  
  pManager = G4SigmaPlus::SigmaPlus()->GetProcessManager();
  
  G4VContinuousDiscreteProcess* theSigmaPlusMult = m_msFactory();
  G4hIonisation* theSigmaPlusIonisation = new G4hIonisation();
  pManager->AddProcess(theSigmaPlusIonisation, ordInActive,2, 2);
  pManager->AddProcess(theSigmaPlusMult);
  pManager->SetProcessOrdering(theSigmaPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theSigmaPlusMult, idxPostStep, 1);

  // AntiSigmaPlus
  
  pManager = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();

  G4VContinuousDiscreteProcess* theAntiSigmaPlusMult = m_msFactory();
  G4hIonisation* theAntiSigmaPlusIonisation = new G4hIonisation();
  pManager->AddProcess(theAntiSigmaPlusIonisation, ordInActive,2, 2);
  pManager->AddProcess(theAntiSigmaPlusMult);
  pManager->SetProcessOrdering(theAntiSigmaPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theAntiSigmaPlusMult, idxPostStep, 1);
  
  // XiMinus
  
  pManager = G4XiMinus::XiMinus()->GetProcessManager();
  
  G4VContinuousDiscreteProcess* theXiMinusMult = m_msFactory() ;
  G4hIonisation* theXiMinusIonisation = new G4hIonisation();
  pManager->AddProcess(theXiMinusIonisation, ordInActive,2, 2);
  pManager->AddProcess(theXiMinusMult);
  pManager->SetProcessOrdering(theXiMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theXiMinusMult, idxPostStep, 1);

  // AntiXiMinus
      
  pManager = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();

  G4VContinuousDiscreteProcess* theAntiXiMinusMult = m_msFactory();
  G4hIonisation* theAntiXiMinusIonisation = new G4hIonisation();
  pManager->AddProcess(theAntiXiMinusIonisation, ordInActive,2, 2);
  pManager->AddProcess(theAntiXiMinusMult);
  pManager->SetProcessOrdering(theAntiXiMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theAntiXiMinusMult, idxPostStep, 1); 

 // OmegaMinus
      
  pManager = G4OmegaMinus::OmegaMinus()->GetProcessManager();
  
  G4VContinuousDiscreteProcess* theOmegaMinusMult = m_msFactory();
  G4hIonisation* theOmegaMinusIonisation = new G4hIonisation();
  pManager->AddProcess(theOmegaMinusIonisation, ordInActive,2, 2);
  pManager->AddProcess(theOmegaMinusMult);
  pManager->SetProcessOrdering(theOmegaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theOmegaMinusMult, idxPostStep, 1); 

// AntiOmegaMinus
      
  pManager = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();

  G4VContinuousDiscreteProcess* theAntiOmegaMinusMult = m_msFactory();
  G4hIonisation* theAntiOmegaMinusIonisation = new G4hIonisation();
  pManager->AddProcess(theAntiOmegaMinusIonisation, ordInActive,2, 2);
  pManager->AddProcess(theAntiOmegaMinusMult);
  pManager->SetProcessOrdering(theAntiOmegaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theAntiOmegaMinusMult, idxPostStep, 1);
  
  // Full Physics List
  
  if (m_physicsChoice=="full" || m_physicsChoice=="improved")
    {
      // Elastic scatter 
      
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess();
      G4LElastic*            theElasticModel;
      
      theElasticModel = new G4LElastic();
      theElasticProcess->RegisterMe(theElasticModel);
  
      // PionPlus
      
      pManager = G4PionPlus::PionPlus()->GetProcessManager();
      
      pManager->AddDiscreteProcess(theElasticProcess);
      G4PionPlusInelasticProcess* thePionPlusInelastic = new G4PionPlusInelasticProcess();
      G4LEPionPlusInelastic* theLEPionPlusModel;
      G4HEPionPlusInelastic* theHEPionPlusModel;
      theLEPionPlusModel = new G4LEPionPlusInelastic();
      theHEPionPlusModel = new G4HEPionPlusInelastic();
      thePionPlusInelastic->RegisterMe(theLEPionPlusModel);
      thePionPlusInelastic->RegisterMe(theHEPionPlusModel);
      pManager->AddDiscreteProcess(thePionPlusInelastic);
    
      // PionMinus
      
      pManager = G4PionMinus::PionMinus()->GetProcessManager();
 
      pManager->AddDiscreteProcess(theElasticProcess);

      G4PionMinusInelasticProcess* thePionMinusInelastic = new G4PionMinusInelasticProcess();
      G4LEPionMinusInelastic* theLEPionMinusModel;
      G4HEPionMinusInelastic* theHEPionMinusModel;
      theLEPionMinusModel = new G4LEPionMinusInelastic();
      theHEPionMinusModel = new G4HEPionMinusInelastic();
      thePionMinusInelastic->RegisterMe(theLEPionMinusModel);
      thePionMinusInelastic->RegisterMe(theHEPionMinusModel);
      pManager->AddDiscreteProcess(thePionMinusInelastic);
      /*
	#ifdef TRIUMF_STOP_PIMINUS
	G4PionMinusAbsorptionAtRest thePionMinusAbsorption;
	#else
	G4PiMinusAbsorptionAtRest thePionMinusAbsorption;
	#endif
	pManager->AddRestProcess(thePionMinusAbsorption, ordDefault);
      */

      // KaonPlus
      
      pManager = G4KaonPlus::KaonPlus()->GetProcessManager();

      pManager->AddDiscreteProcess(theElasticProcess);
      G4KaonPlusInelasticProcess* theKaonPlusInelastic = new G4KaonPlusInelasticProcess();
      G4LEKaonPlusInelastic* theLEKaonPlusModel;
      G4HEKaonPlusInelastic* theHEKaonPlusModel;      
      theLEKaonPlusModel = new G4LEKaonPlusInelastic();
      theHEKaonPlusModel = new G4HEKaonPlusInelastic();
      theKaonPlusInelastic->RegisterMe(theLEKaonPlusModel);
      theKaonPlusInelastic->RegisterMe(theHEKaonPlusModel);
      pManager->AddDiscreteProcess(theKaonPlusInelastic);
  
      // KaonMinus
      
      pManager = G4KaonMinus::KaonMinus()->GetProcessManager();
      
      pManager->AddDiscreteProcess(theElasticProcess);
      G4KaonMinusInelasticProcess* theKaonMinusInelastic = new G4KaonMinusInelasticProcess();
      G4LEKaonMinusInelastic* theLEKaonMinusModel;
      G4HEKaonMinusInelastic* theHEKaonMinusModel;
      theLEKaonMinusModel = new G4LEKaonMinusInelastic();
      theHEKaonMinusModel = new G4HEKaonMinusInelastic();
      theKaonMinusInelastic->RegisterMe(theLEKaonMinusModel);
      theKaonMinusInelastic->RegisterMe(theHEKaonMinusModel);
      pManager->AddDiscreteProcess(theKaonMinusInelastic);

      /*
	#ifdef TRIUMF_STOP_KMINUS
	G4KaonMinusAbsorption* theKaonMinusAbsorption = new G4KaonMinusAbsorption();
	#else
	G4PiMinusAbsorptionAtRest* theKaonMinusAbsorption = new G4PiMinusAbsorptionAtRest();
	#endif
	pManager->AddRestProcess(theKaonMinusAbsorption, ordDefault);
      */

      // KaonZeroL
      
      pManager = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
      
      pManager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroLInelasticProcess* theKaonZeroLInelastic = new G4KaonZeroLInelasticProcess();
      G4LEKaonZeroLInelastic* theLEKaonZeroLModel;
      G4HEKaonZeroInelastic* theHEKaonZeroLModel;
      theLEKaonZeroLModel = new G4LEKaonZeroLInelastic();
      theHEKaonZeroLModel = new G4HEKaonZeroInelastic();
      theKaonZeroLInelastic->RegisterMe(theLEKaonZeroLModel);
      theKaonZeroLInelastic->RegisterMe(theHEKaonZeroLModel);
      pManager->AddDiscreteProcess(theKaonZeroLInelastic);
  
      // KaonZeroS
  
      pManager = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
  
      pManager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroSInelasticProcess* theKaonZeroSInelastic = new G4KaonZeroSInelasticProcess(); 
      G4LEKaonZeroSInelastic* theLEKaonZeroSModel;
      G4HEKaonZeroInelastic* theHEKaonZeroSModel;
      theLEKaonZeroSModel = new G4LEKaonZeroSInelastic();
      theHEKaonZeroSModel = new G4HEKaonZeroInelastic();
      theKaonZeroSInelastic->RegisterMe(theLEKaonZeroSModel);
      theKaonZeroSInelastic->RegisterMe(theHEKaonZeroSModel);
      pManager->AddDiscreteProcess(theKaonZeroSInelastic);

      // Proton
      
      pManager = G4Proton::Proton()->GetProcessManager();

      pManager->AddDiscreteProcess(theElasticProcess);  
      G4ProtonInelasticProcess* theProtonInelastic = new G4ProtonInelasticProcess();  
      G4LEProtonInelastic* theLEProtonModel;
      G4HEProtonInelastic* theHEProtonModel;
      theLEProtonModel = new G4LEProtonInelastic();
      theHEProtonModel = new G4HEProtonInelastic();
      theProtonInelastic->RegisterMe(theLEProtonModel);
      theProtonInelastic->RegisterMe(theHEProtonModel);
      pManager->AddDiscreteProcess(theProtonInelastic);

      // AntiProton
  
      pManager = G4AntiProton::AntiProton()->GetProcessManager();

      pManager->AddDiscreteProcess(theElasticProcess);
      G4AntiProtonInelasticProcess* theAntiProtonInelastic = new G4AntiProtonInelasticProcess(); 
      G4LEAntiProtonInelastic* theLEAntiProtonModel;
      G4HEAntiProtonInelastic* theHEAntiProtonModel;   
      theLEAntiProtonModel = new G4LEAntiProtonInelastic();
      theHEAntiProtonModel = new G4HEAntiProtonInelastic();
      theAntiProtonInelastic->RegisterMe(theLEAntiProtonModel);
      theAntiProtonInelastic->RegisterMe(theHEAntiProtonModel);
      pManager->AddDiscreteProcess(theAntiProtonInelastic);

      // Neutron

      pManager = G4Neutron::Neutron()->GetProcessManager();

      pManager->AddDiscreteProcess(theElasticProcess);
      G4NeutronInelasticProcess*  theNeutronInelastic = new G4NeutronInelasticProcess();
      G4LENeutronInelastic* theLENeutronModel;
      G4HENeutronInelastic* theHENeutronModel;
      theLENeutronModel = new G4LENeutronInelastic();
      theHENeutronModel = new G4HENeutronInelastic();
      theNeutronInelastic->RegisterMe(theLENeutronModel);
      theNeutronInelastic->RegisterMe(theHENeutronModel);
      pManager->AddDiscreteProcess(theNeutronInelastic);
      /* 
	 G4HadronFissionProcess* theNeutronFission = new G4HadronFissionProcess();
	 G4LFission* theNeutronFissionModel;
	 G4HadronCaptureProcess* theNeutronCapture = new G4HadronCaptureProcess();
	 G4LCapture* theNeutronCaptureModel;
	 theNeutronFissionModel = new G4LFission();
	 theNeutronFission->RegisterMe(theNeutronFissionModel);
	 pManager->AddDiscreteProcess(theNeutronFission);
	 theNeutronCaptureModel = new G4LCapture();
	 theNeutronCapture->RegisterMe(theNeutronCaptureModel);
	 pManager->AddDiscreteProcess(theNeutronCapture);
      */

      // AntiNeutron

      pManager = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  
      pManager->AddDiscreteProcess(theElasticProcess);
      G4AntiNeutronInelasticProcess*  theAntiNeutronInelastic = new G4AntiNeutronInelasticProcess();
      G4LEAntiNeutronInelastic* theLEAntiNeutronModel;
      G4HEAntiNeutronInelastic* theHEAntiNeutronModel;      
      theLEAntiNeutronModel = new G4LEAntiNeutronInelastic();
      theHEAntiNeutronModel = new G4HEAntiNeutronInelastic();
      theAntiNeutronInelastic->RegisterMe(theLEAntiNeutronModel);
      theAntiNeutronInelastic->RegisterMe(theHEAntiNeutronModel);
      pManager->AddDiscreteProcess(theAntiNeutronInelastic);
     
      // Lambda
      
      pManager = G4Lambda::Lambda()->GetProcessManager();
      
      pManager->AddDiscreteProcess(theElasticProcess);
      G4LambdaInelasticProcess*  theLambdaInelastic = new G4LambdaInelasticProcess();
      G4LELambdaInelastic*  theLELambdaModel;
      G4HELambdaInelastic*  theHELambdaModel;
      theLELambdaModel = new G4LELambdaInelastic();
      theHELambdaModel = new G4HELambdaInelastic();
      theLambdaInelastic->RegisterMe(theLELambdaModel);
      theLambdaInelastic->RegisterMe(theHELambdaModel);
      pManager->AddDiscreteProcess(theLambdaInelastic);
  
      // AntiLambda
      
      pManager = G4AntiLambda::AntiLambda()->GetProcessManager();
  
      pManager->AddDiscreteProcess(theElasticProcess);
      G4AntiLambdaInelasticProcess*  theAntiLambdaInelastic = new G4AntiLambdaInelasticProcess();
      G4LEAntiLambdaInelastic*  theLEAntiLambdaModel;
      G4HEAntiLambdaInelastic*  theHEAntiLambdaModel;
      theLEAntiLambdaModel = new G4LEAntiLambdaInelastic();
      theHEAntiLambdaModel = new G4HEAntiLambdaInelastic();
      theAntiLambdaInelastic->RegisterMe(theLEAntiLambdaModel);
      theAntiLambdaInelastic->RegisterMe(theHEAntiLambdaModel);
      pManager->AddDiscreteProcess(theAntiLambdaInelastic);
      
      // SigmaMinus
      
      pManager = G4SigmaMinus::SigmaMinus()->GetProcessManager();

      pManager->AddDiscreteProcess(theElasticProcess);
      G4SigmaMinusInelasticProcess*  theSigmaMinusInelastic = new G4SigmaMinusInelasticProcess(); 
      G4LESigmaMinusInelastic*  theLESigmaMinusModel;
      G4HESigmaMinusInelastic*  theHESigmaMinusModel;
      theLESigmaMinusModel = new G4LESigmaMinusInelastic();
      theHESigmaMinusModel = new G4HESigmaMinusInelastic();
      theSigmaMinusInelastic->RegisterMe(theLESigmaMinusModel);
      theSigmaMinusInelastic->RegisterMe(theHESigmaMinusModel);
      pManager->AddDiscreteProcess(theSigmaMinusInelastic);

      // AntiSigmaMinus
      
      pManager = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
      
      pManager->AddDiscreteProcess(theElasticProcess);
      G4AntiSigmaMinusInelasticProcess*  theAntiSigmaMinusInelastic = new G4AntiSigmaMinusInelasticProcess();
      G4LEAntiSigmaMinusInelastic*  theLEAntiSigmaMinusModel;
      G4HEAntiSigmaMinusInelastic*  theHEAntiSigmaMinusModel;
      theLEAntiSigmaMinusModel = new G4LEAntiSigmaMinusInelastic();
      theHEAntiSigmaMinusModel = new G4HEAntiSigmaMinusInelastic();
      theAntiSigmaMinusInelastic->RegisterMe(theLEAntiSigmaMinusModel);
      theAntiSigmaMinusInelastic->RegisterMe(theHEAntiSigmaMinusModel);
      pManager->AddDiscreteProcess(theAntiSigmaMinusInelastic);
 
      // SigmaPlus
      
      pManager = G4SigmaPlus::SigmaPlus()->GetProcessManager();
      
      pManager->AddDiscreteProcess(theElasticProcess);
      G4SigmaPlusInelasticProcess*  theSigmaPlusInelastic = new G4SigmaPlusInelasticProcess();
      G4LESigmaPlusInelastic*  theLESigmaPlusModel;
      G4HESigmaPlusInelastic*  theHESigmaPlusModel;  
      theLESigmaPlusModel = new G4LESigmaPlusInelastic();
      theHESigmaPlusModel = new G4HESigmaPlusInelastic();
      theSigmaPlusInelastic->RegisterMe(theLESigmaPlusModel);
      theSigmaPlusInelastic->RegisterMe(theHESigmaPlusModel);
      pManager->AddDiscreteProcess(theSigmaPlusInelastic);

      // AntiSigmaPlus
      
      pManager = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();

      pManager->AddDiscreteProcess(theElasticProcess);
      G4AntiSigmaPlusInelasticProcess*  theAntiSigmaPlusInelastic = new G4AntiSigmaPlusInelasticProcess();
      G4LEAntiSigmaPlusInelastic*  theLEAntiSigmaPlusModel;
      G4HEAntiSigmaPlusInelastic*  theHEAntiSigmaPlusModel;   
      theLEAntiSigmaPlusModel = new G4LEAntiSigmaPlusInelastic();
      theHEAntiSigmaPlusModel = new G4HEAntiSigmaPlusInelastic();
      theAntiSigmaPlusInelastic->RegisterMe(theLEAntiSigmaPlusModel);
      theAntiSigmaPlusInelastic->RegisterMe(theHEAntiSigmaPlusModel);
      pManager->AddDiscreteProcess(theAntiSigmaPlusInelastic);
    
      // XiZero
      
      pManager = G4XiZero::XiZero()->GetProcessManager();
      
      pManager->AddDiscreteProcess(theElasticProcess);
      G4XiZeroInelasticProcess*  theXiZeroInelastic = new G4XiZeroInelasticProcess();
      G4LEXiZeroInelastic*  theLEXiZeroModel;
      G4HEXiZeroInelastic*  theHEXiZeroModel;    
      theLEXiZeroModel = new G4LEXiZeroInelastic();
      theHEXiZeroModel = new G4HEXiZeroInelastic();
      theXiZeroInelastic->RegisterMe(theLEXiZeroModel);
      theXiZeroInelastic->RegisterMe(theHEXiZeroModel);
      pManager->AddDiscreteProcess(theXiZeroInelastic);

      // AntiXiZero

      pManager = G4AntiXiZero::AntiXiZero()->GetProcessManager();

      pManager->AddDiscreteProcess(theElasticProcess);
      G4AntiXiZeroInelasticProcess*  theAntiXiZeroInelastic = new G4AntiXiZeroInelasticProcess();
      G4LEAntiXiZeroInelastic*  theLEAntiXiZeroModel;
      G4HEAntiXiZeroInelastic*  theHEAntiXiZeroModel;
      theLEAntiXiZeroModel = new G4LEAntiXiZeroInelastic();
      theHEAntiXiZeroModel = new G4HEAntiXiZeroInelastic();
      theAntiXiZeroInelastic->RegisterMe(theLEAntiXiZeroModel);
      theAntiXiZeroInelastic->RegisterMe(theHEAntiXiZeroModel);
      pManager->AddDiscreteProcess(theAntiXiZeroInelastic);
  
      // XiMinus
  
      pManager = G4XiMinus::XiMinus()->GetProcessManager();
  
      pManager->AddDiscreteProcess(theElasticProcess);
      G4XiMinusInelasticProcess*  theXiMinusInelastic = new G4XiMinusInelasticProcess();
      G4LEXiMinusInelastic*  theLEXiMinusModel;
      G4HEXiMinusInelastic*  theHEXiMinusModel;
      theLEXiMinusModel = new G4LEXiMinusInelastic();
      theHEXiMinusModel = new G4HEXiMinusInelastic();
      theXiMinusInelastic->RegisterMe(theLEXiMinusModel);
      theXiMinusInelastic->RegisterMe(theHEXiMinusModel);
      pManager->AddDiscreteProcess(theXiMinusInelastic);

      // AntiXiMinus
      
      pManager = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();

      pManager->AddDiscreteProcess(theElasticProcess);
      G4AntiXiMinusInelasticProcess*  theAntiXiMinusInelastic = new G4AntiXiMinusInelasticProcess();
      G4LEAntiXiMinusInelastic*  theLEAntiXiMinusModel;
      G4HEAntiXiMinusInelastic*  theHEAntiXiMinusModel;  
      theLEAntiXiMinusModel = new G4LEAntiXiMinusInelastic();
      theHEAntiXiMinusModel = new G4HEAntiXiMinusInelastic();
      theAntiXiMinusInelastic->RegisterMe(theLEAntiXiMinusModel);
      theAntiXiMinusInelastic->RegisterMe(theHEAntiXiMinusModel);
      pManager->AddDiscreteProcess(theAntiXiMinusInelastic);
 
      // OmegaMinus
      
      pManager = G4OmegaMinus::OmegaMinus()->GetProcessManager();
      
      pManager->AddDiscreteProcess(theElasticProcess);
      G4OmegaMinusInelasticProcess*  theOmegaMinusInelastic = new G4OmegaMinusInelasticProcess();
      G4LEOmegaMinusInelastic*  theLEOmegaMinusModel;
      G4HEOmegaMinusInelastic*  theHEOmegaMinusModel;
      theLEOmegaMinusModel = new G4LEOmegaMinusInelastic();
      theHEOmegaMinusModel = new G4HEOmegaMinusInelastic();
      theOmegaMinusInelastic->RegisterMe(theLEOmegaMinusModel);
      theOmegaMinusInelastic->RegisterMe(theHEOmegaMinusModel);
      pManager->AddDiscreteProcess(theOmegaMinusInelastic);

      // AntiOmegaMinus
      
      pManager = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();

      pManager->AddDiscreteProcess(theElasticProcess);
      G4AntiOmegaMinusInelasticProcess*  theAntiOmegaMinusInelastic = new G4AntiOmegaMinusInelasticProcess();
      G4LEAntiOmegaMinusInelastic*  theLEAntiOmegaMinusModel;
      G4HEAntiOmegaMinusInelastic*  theHEAntiOmegaMinusModel; 
      theLEAntiOmegaMinusModel = new G4LEAntiOmegaMinusInelastic();
      theHEAntiOmegaMinusModel = new G4HEAntiOmegaMinusInelastic();
      theAntiOmegaMinusInelastic->RegisterMe(theLEAntiOmegaMinusModel);
      theAntiOmegaMinusInelastic->RegisterMe(theHEAntiOmegaMinusModel);
      pManager->AddDiscreteProcess(theAntiOmegaMinusInelastic);
    
    }
  
}












