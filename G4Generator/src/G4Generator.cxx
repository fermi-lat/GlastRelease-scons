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

#include "RunManager.h"
#include "PrimaryGeneratorAction.h"
#include "McParticleManager.h"
// GLAST Geant4

#include "GlastMS/MultipleScatteringFactory.h"  // For controlling multiple scattering versions
#include "GlastMS/EnergyLossFactory.h"          // For controlling energy loss versions in EMPhysics

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
#include "Event/MonteCarlo/McTrajectory.h"

#include "G4Generator/IG4GeometrySvc.h"

//vectors
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"


static const AlgFactory<G4Generator>  Factory;
const IAlgFactory& G4GeneratorFactory = Factory;

namespace {
    double zOffset;
}

G4Generator::G4Generator(const std::string& name, ISvcLocator* pSvcLocator) 
  :Algorithm(name, pSvcLocator) 
{
  // set defined properties
  declareProperty("UIcommands", m_uiCommands);
  declareProperty("geometryMode", m_geometryMode="");
  declareProperty("saveTrajectories", m_saveTrajectories=0);
  declareProperty("mcTreeMode", m_mcTreeMode="minimal");
  declareProperty("defaultCutValue", m_defaultCutValue=0.1*mm);
  declareProperty("defaultTkrCutValue", m_defaultTkrCutValue=0.1*mm);
  declareProperty("defaultCalCutValue", m_defaultCalCutValue=0.1*mm);
  declareProperty("physics_choice", m_physics_choice="full");
  declareProperty("physics_tables", m_physics_table="build");
  declareProperty("physics_dir", m_physics_dir="G4cuts/100micron/");
  declareProperty("numGenerations", m_numGenerations=100000);
  declareProperty("lowEnergyCut",m_lowEnergy=0.52);
  declareProperty("minDistanceCut",m_minDistance=0.1);

  declareProperty("mscatOption",  m_mscatOption  = true); 
  declareProperty("eLossCurrent", m_eLossCurrent = true);  // default is to use the current, not 5.2.
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

  if( service( "ParticlePropertySvc", m_ppsvc).isFailure() ) {
      log << MSG::ERROR << "Couldn't set up ParticlePropertySvc!" << endreq;
      return StatusCode::FAILURE;
  }

  // create a factory for the multiplescattering to pass around to the physics guys
  GlastMS::MultipleScatteringFactory msFactory(
      m_mscatOption ? GlastMS::MultipleScatteringFactory::OLD32 
                    : GlastMS::MultipleScatteringFactory::NATIVE);
  GlastMS::EnergyLossFactory eLossFactory(
      m_eLossCurrent ? GlastMS::EnergyLossFactory::CURRENT 
                     : GlastMS::EnergyLossFactory::RELEASE52);

  log << MSG::INFO << "Using the " << (m_mscatOption? "Old 3.2" : "current G4") 
      << " version of Multiple scattering" << endreq;

  log << MSG::INFO << "Using the " << (!m_eLossCurrent? "Old 5.2" : "current G4") 
      << " version of Energy Loss" << endreq;
  
  log << MSG::INFO << "Initializing run manager...\n";
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
                                    msFactory,
                                    eLossFactory,
									geosv);

      log << "\n done." << endreq;
    }

  // Initialize Geant4
  m_runManager->Initialize();

  // set the mode for the McParticle tree
  if (m_mcTreeMode == "minimal")
    McParticleManager::getPointer()->setMode(0);
 
  log << endreq;
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

  // Clean the McParticle hierarchy 
  McParticleManager::getPointer()->clear();

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

  assert(pcol->size()==1); // something wrong: must be only one
  Event::McParticle* primary = pcol->front();

  // we get the primaryGenerator from the RunManager
  PrimaryGeneratorAction* primaryGenerator = 
    (PrimaryGeneratorAction*)m_runManager->GetUserPrimaryGeneratorAction();
  
  // we set the primary particle by passing the pointer of the McParticle to the
  // primaryGenerator class
  primaryGenerator->init(primary, m_ppsvc, zOffset);

  McParticleManager::getPointer()->addMcParticle(0,primary);   
  // Run geant4
  m_runManager->BeamOn(); 
    
  // If we set the pruneCal mode than we remove all the not TKR interacting
  // particles that does not touch the TKR zone (z>0)
  if (m_mcTreeMode=="pruneCal")
    McParticleManager::getPointer()->pruneCal();

  // Save the McParticle hierarchy in the TDS
  McParticleManager::getPointer()->save();

  // we save the trajectories if the m_saveTrajectories is true
  if (m_saveTrajectories)
    {
      Event::McTrajectoryCol* traj = new Event::McTrajectoryCol();  
      eventSvc()->registerObject("/Event/MC/TrajectoryCol",traj);

      int numTrajs  = m_runManager->getNumberOfTrajectories();
      int maxGenNum = 0;

      std::multimap<const int, const unsigned int> genVsTrajIdx;

      for(unsigned int j = 0; j < m_runManager->getNumberOfTrajectories(); j++)
      {
          // Start by finding the McParticle associated with this trajectory
          // note that even if we can't find an associated particle, we store the trajectory and its sign
          Event::McParticle* part = 0;

          if (McParticleManager::getPointer()->size() > 
              static_cast<unsigned int>( m_runManager->getTrajectoryTrackId(j)))
          {
              part = 
                  McParticleManager::getPointer()->getMcParticle(m_runManager->getTrajectoryTrackId(j));
          }

          int genNum = 1000;

          if (part)
          {
              if (part->initialFourMomentum().e() < m_lowEnergy)
              {
                  double dist = (part->initialPosition() - part->finalPosition()).mag();

                  if (dist < m_minDistance) continue;
              }

              genNum = 0;
              SmartRef<Event::McParticle> mcPart = part->getMother();

              while(mcPart != mcPart->getMother())
              {
                  mcPart = mcPart->getMother();
                  genNum++;
              }

              if (genNum > maxGenNum) maxGenNum = genNum;
          }

          genVsTrajIdx.insert(std::pair<const int, const unsigned int>(genNum,j));
      }

      // Check size of map
      int genIdxMapSize  = genVsTrajIdx.size();
      int maxGenerations = 10000;

      if (genIdxMapSize > 20000)
      {
          maxGenerations = maxGenNum - 2;
          if (maxGenerations > m_numGenerations) maxGenerations = m_numGenerations;
      }

      // Now go through the map and set up generation related display
      std::multimap<const int, const unsigned int>::iterator genMapIter;
      for(genMapIter = genVsTrajIdx.begin(); genMapIter != genVsTrajIdx.end(); genMapIter++)
      {
          const int          genNum  = (*genMapIter).first;
          const unsigned int trajIdx = (*genMapIter).second;

          if (genNum <= maxGenerations)
          {
              Event::McParticle* part = 0;

              if (McParticleManager::getPointer()->size() > 
                  static_cast<unsigned int>( m_runManager->getTrajectoryTrackId(trajIdx)))
              {
                  part = 
                      McParticleManager::getPointer()->getMcParticle(m_runManager->getTrajectoryTrackId(trajIdx));
              }

              std::auto_ptr<std::vector<Hep3Vector> > points = 
                  m_runManager->getTrajectoryPoints(trajIdx);
              int numPoints = (*(points.get())).size();
              Event::McTrajectory* tr = new Event::McTrajectory();
              tr->addPoints(*(points.get()));
              tr->setCharge(m_runManager->getTrajectoryCharge(trajIdx));
              tr->setMcParticle(part); //may be zero.

              traj->push_back(tr);
          }
      }
  }

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





