/** @file G4Generator.cxx
    @brief implementation of class G4Generator

    $Header$

 This is the Gaudi algorithm that runs Geant4 and fills the TDS
 with Montecarlo data. It initalizes some services (for tds and detector
 geometry) and than passes them to the RunManager class that do most of the
 work. Optionally it can uses the GuiSvc to show the event in a GUI; in that
 case the GUI takes control on the event loop
      
 @author  T.Burnett
 @author  R.Giannitrapani
*/

// Include files

// Geant4
#include "G4Generator.h"
#include "G4UImanager.hh"
#include "G4ParticleDefinition.hh"

#include "G4Material.hh"

#include "RunManager.h"
#include "PrimaryGeneratorAction.h"
#include "McParticleManager.h"
#include "McTrajectoryManager.h"

// Gaudi
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/IAlgManager.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

//special to setup the TdGlastData structure
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

//montecarlo data structures 
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
//#include "Event/MonteCarlo/McTrajectory.h"

#include "G4Generator/IG4GeometrySvc.h"

//vectors
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"

// Error stuff
#include "G4Generator/IG4GenErrorSvc.h"
#include "G4Generator/G4GenException.h"

// General
#include <iomanip>

static const AlgFactory<G4Generator>  Factory;
const IAlgFactory& G4GeneratorFactory = Factory;

namespace {
    double zOffset;
}

G4Generator::G4Generator(const std::string& name, ISvcLocator* pSvcLocator) 
  :Algorithm(name, pSvcLocator) 
{
  // set defined properties
  declareProperty("UIcommands",   m_uiCommands);
  declareProperty("geometryMode", m_geometryMode="");

  // This is no longer used but kept to preserve old jo files
  declareProperty("saveTrajectories", m_saveTrajectories=0);

  // G4 cut offs for tracking
  declareProperty("defaultCutValue",    m_defaultCutValue=0.1*mm);
  declareProperty("defaultTkrCutValue", m_defaultTkrCutValue=0.1*mm);
  declareProperty("defaultCalCutValue", m_defaultCalCutValue=0.1*mm);

  // G4 choices 
  declareProperty("physics_choice", m_physics_choice="QGSP_BERT");
  declareProperty("physics_tables", m_physics_table="build");
  declareProperty("physics_dir",    m_physics_dir="G4cuts/100micron/");

  // McParticle/McTrajectory pruning modes
  // "full" means no pruning of particles/trajectories
  // "minimal" means primary and direct daughters only (but > m_lowEnergy)
  // "nGenerations" means keep up to m_numGenerations generations of daughters
  // "prunecal" means prune particles originating and termination in the cal
  declareProperty("mcTreeMode",     m_mcTreeMode="minimal");
  declareProperty("numGenerations", m_numGenerations=4);
  declareProperty("lowEnergyCut",   m_lowEnergy=0.52);
  declareProperty("minDistanceCut", m_minDistance=0.1);

  declareProperty("mscatOption",  m_mscatOption  = true); 
  declareProperty("eLossCurrent", m_eLossCurrent = true);  // default is to use the current, not 5.2.

  declareProperty("printRadLen",  m_printRadLen  = true);
  // if printRadLen is true, print all elements, including duplicates
  declareProperty("printAll",     m_printAll     = false);
}
    
////////////////////////////////////////////////////////////////////////////
StatusCode G4Generator::initialize()
{
  // Purpose and Method: This routine initialize the Gaudi algorithm.  It is
  // called once before event processing begins.  
  // Outputs: A StatusCode which denotes success or failure.

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "initialize" << endreq;

  // Use the Job options service to set the Algorithm's parameters
  setProperties();
  // user jobOtions on cuts and physics 

  log << MSG::INFO << "DefaultCutValue="    << m_defaultCutValue    << " mm" << endreq;
  log << MSG::INFO << "DefaultTkrCutValue=" << m_defaultTkrCutValue << " mm" << endreq;
  log << MSG::INFO << "DefaultCalCutValue=" << m_defaultCalCutValue << " mm" << endreq;

  log << MSG::INFO << "Physics List = " << m_physics_choice << endreq;

  // user jobOtions on cuts directory 

  log << MSG::INFO << "Cut Table " << m_physics_table << endreq;
  log << MSG::INFO << "Cut Table Dir" << m_physics_dir << endreq;
 
  // Apply Geant4 specific commands throught the ui
  if( !m_uiCommands.value().empty() ) {
    G4UImanager* UI = G4UImanager::GetUIpointer();
    for( std::vector<std::string>::const_iterator k = 
           m_uiCommands.value().begin(); 
         k!=m_uiCommands.value().end(); ++k){
      UI->ApplyCommand(*k);
      log << MSG::INFO << "UI command: " << (*k) << endreq;
    }
  }  

  // Get the G4Geometry Service
  IG4GeometrySvc* geosv=0;
  if( service( "G4GeometrySvc", geosv, true).isFailure() ) {
    log << MSG::ERROR << "Couldn't set up G4GeometrySvc!" << endreq;
    return StatusCode::FAILURE;
  }

  // Error handling service 
  if (service("G4GenErrorSvc",m_ErrorSvc, true).isFailure())
  {
      log << MSG::ERROR << "Could not find G4GenErrorSvc"<<endreq ;
      return StatusCode::FAILURE ;
  }

  // Get the GlastDetService Service
  IGlastDetSvc* glastsvc=0;
  if( service( "GlastDetSvc", glastsvc, true).isFailure() ) {
    log << MSG::ERROR << "Couldn't set up GlastDetSvc!" << endreq;
    return StatusCode::FAILURE;
  }

  //look for z offset
  zOffset = 0.0;
  if (glastsvc->getTopVolume()=="LAT") {
      if(glastsvc->getNumericConstByName("globalOffset_dz", &zOffset).isFailure()) 
          zOffset = 0.0;
  }
  if (zOffset!=0.0) {
      log << MSG::INFO << "Particles will be offset in z by " 
          << zOffset << " mm" << endreq;
  }

  // Init the McParticle hierarchy 
  McParticleManager::getPointer()->initialize(eventSvc(), glastsvc);

  // Init the McTrajectory hierarchy 
  McTrajectoryManager::getPointer()->initialize(eventSvc(), glastsvc);

  if( service( "ParticlePropertySvc", m_ppsvc).isFailure() ) {
      log << MSG::ERROR << "Couldn't set up ParticlePropertySvc!" << endreq;
      return StatusCode::FAILURE;
  }

  // create a factory for the multiplescattering to pass around to the physics guys
  
  log << MSG::INFO << "Initializing run manager... ";

  // The geant4 manager
  if (!(m_runManager = RunManager::GetRunManager()))
    {
      m_runManager = new RunManager(log.stream(), 
                                    m_defaultCutValue, 
                                    m_defaultTkrCutValue,
                                    m_defaultCalCutValue,
                                    m_physics_choice, 
                                    m_physics_table, 
                                    m_physics_dir,
				    geosv);

      log << "done." << endreq;
    }

  // Initialize Geant4
  m_runManager->Initialize();

  // set the mode for the McParticle tree
  if (m_mcTreeMode == "nGenerations")
  {
      McParticleManager::getPointer()->setMode(McParticleManager::N_GENERATIONS);
      McParticleManager::getPointer()->setNumGenerations(m_numGenerations);
  }
  else if (m_mcTreeMode == "prunecal")
  {
      McParticleManager::getPointer()->setMode(McParticleManager::PRUNE_CAL);
  }
  else if (m_mcTreeMode == "full")
  {
      McParticleManager::getPointer()->setMode(McParticleManager::FULL_MC_TREE);
  }
  else // default mode is "minimal" 
  {
      McParticleManager::getPointer()->setMode(McParticleManager::MINIMAL_TREE);
      McParticleManager::getPointer()->setCutOffEnergy(m_lowEnergy);
  }
  log << MSG::INFO << endreq  
      << "List of materials and radiation lengths" 
      << endreq;

  if(m_printRadLen) {
      std::cout << std::endl << std::left
          << std::setw(10) << "Index"
          << std::setw(16) << "Name" 
          << std::setw(9) << "Density" 
          << std::setw(15) << "RadLength (cm)" 
          << std::setw(8) << "(g/cm3)" << std::endl << std::endl;

      int size = G4Material::GetNumberOfMaterials(); 
      int i;
      for(i=0;i<size;++i) {
          G4Material* mat = (*G4Material::GetMaterialTable())[i];
          double density = mat->GetDensity()/(g/cm3);
          double radlen  = mat->GetRadlen()/cm;
          double radlenG = radlen*density;
          G4Material* matName = G4Material::GetMaterial(mat->GetName());
          int indName = matName->GetIndex();
          if(!m_printAll&&i!=indName) continue;
          std::cout << std::left << std::setprecision(4) 
              << std::setw(2) << i; 
          //if(i!=indName) {
          //    std::cout << "(" << std::setw(2) << indName << ") ";
          //} else {
              std::cout << "     "; 
          //}
          std::cout << std::left   << std::setw(20) << mat->GetName() 
              << std::setw(12) << density
              << std::setw(13) << radlen
              << std::setw(8) << radlenG << std::endl;
      } 
      std::cout << std::endl;
  } else {
      log << "... Not requested" << endreq;
  }
      //std::cout << "\n\n\n";
      //std::cout << *(G4Material::GetMaterialTable()) << std::endl;

  return StatusCode::SUCCESS;
}


StatusCode G4Generator::execute() 
{
  // Purpose and Method: This is the execute routine of the algorithm, called
  //            once every event. 
  // Outputs:  A StatusCode which denotes success or failure.
  // TDS Inputs:  If it exists, the McParticle tree to get the primary

  MsgStream   log( msgSvc(), name() );

  // Here the TDS is prepared to receive hits vectors
  // Check for the MC branch - it will be created if it is not available
  DataObject *mc;
  eventSvc()->retrieveObject("/Event/MC", mc);

  log << MSG::DEBUG << "TDS ready" << endreq;

  // Clean the McParticle & McTrajectory hierarchy 
  McParticleManager::getPointer()->clear();
  McTrajectoryManager::getPointer()->clear();

  // following model from previous version, allow property "UIcommands" to
  // generate UI commands here.
  //
  if( !m_uiCommands.value().empty() ) {
    for( std::vector<std::string>::const_iterator 
           k = m_uiCommands.value().begin(); 
         k!=m_uiCommands.value().end(); ++k){
      G4UImanager::GetUIpointer()->ApplyCommand(*k);
      log << MSG::DEBUG << "Apply UI command: \"" << (*k) << "\"" <<endreq;
    }
  }  

  Event::McParticleCol*  pcol=  
    SmartDataPtr<Event::McParticleCol>(eventSvc(), "/Event/MC/McParticleCol");

  if( pcol==0){ 
      log<< MSG::ERROR << "No source of particles!" << endreq;
      return StatusCode::FAILURE;
  
  }  

  assert(pcol->size() > 0); // something wrong: must be 1 or more

  // we get the primaryGenerator from the RunManager
  PrimaryGeneratorAction* primaryGenerator = 
    (PrimaryGeneratorAction*)m_runManager->GetUserPrimaryGeneratorAction();

  // Put the following in a try-catch block to catch errors which might happen during this phase
  try 
  {
      // we set the primary particle by passing the pointer of the McParticle to the
      // primaryGenerator class
      primaryGenerator->init(pcol, m_ppsvc, zOffset);

      // Add the primary particles to the McParticle list
      McParticleManager::getPointer()->addMcParticle(pcol);

      // Run geant4
      m_runManager->BeamOn(); 
  }
  catch(G4GenException& e)
  {
      StatusCode sc = m_ErrorSvc->handleError(name()+" G4Exception",e.what());
      if (sc == StatusCode::SUCCESS) m_runManager->RunTermination();
      return sc;
  }
  catch(...)
  {
      StatusCode sc = m_ErrorSvc->handleError(name(),"unknown exception");
      if (sc == StatusCode::SUCCESS) m_runManager->RunTermination();
      return sc;
  }
    
  // If we set the pruneCal mode than we remove all the not TKR interacting
  // particles that does not touch the TKR zone (z>0)
  if (m_mcTreeMode=="pruneCal")
    McParticleManager::getPointer()->pruneCal();

  // Save the McParticle hierarchy in the TDS
  McParticleManager::getPointer()->save();
  McTrajectoryManager::getPointer()->save();

  return StatusCode::SUCCESS;
}

StatusCode G4Generator::finalize() 
{
  // Purpose and Method: This is the finalize routine of the algorithm, called
  //            at the end of the event loop. 
  // Outputs:  A StatusCode which denotes success or failure.

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "finalize: " << endreq;

  // delete the runManager of geant4
  delete m_runManager;
    
  return StatusCode::SUCCESS;
}





