// File and Version Information:
// $Header$
//
// Description: This is the Gaudi algorithm that runs Geant4 and fills the TDS
// with Montecarlo data. It initalizes some services (for tds and detector
// geometry) and than passes them to the RunManager class that do most of the
// work. Optionally it can uses the GuiSvc to show the event in a GUI; in that
// case the GUI takes control on the event loop
//      
// Author(s):
//      T.Burnett
//      R.Giannitrapani


// Include files

// Geant4
#include "G4Generator.h"
#include "G4UImanager.hh"
#include "G4ParticleDefinition.hh"

#include "RunManager.h"
#include "PrimaryGeneratorAction.h"
#include "McParticleManager.h"

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

//gui
#include "GuiSvc/GuiSvc.h"
#include "gui/DisplayControl.h"
#include "gui/GuiMgr.h"
#include "gui/SimpleCommand.h"
#include "DisplayManager.h"

#include "G4Generator/IG4GeometrySvc.h"

//vectors
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"


static const AlgFactory<G4Generator>  Factory;
const IAlgFactory& G4GeneratorFactory = Factory;

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

  // setup the GuiSvc, if available
  setupGui();

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

  // Init the McParticle hierarchy 
  McParticleManager::getPointer()->initialize(eventSvc());

  if( service( "ParticlePropertySvc", m_ppsvc).isFailure() ) {
      log << MSG::ERROR << "Couldn't set up ParticlePropertySvc!" << endreq;
      return StatusCode::FAILURE;
  }

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
                                    m_physics_dir);

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

void G4Generator::setupGui()
{
  // Purpose and Method:  This routine setup the (optional) Gui service

  MsgStream log(msgSvc(), name());

  IGuiSvc* guiSvc=0;
    
  if ( service("GuiSvc", guiSvc).isFailure() ){
    log << MSG::WARNING << "No GuiSvc: so, no event display " << endreq;
    return;
  } 
  new DisplayManager(&(guiSvc->guiMgr()->display()));
   
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
      log << MSG::INFO << "Apply UI command: \"" << (*k) << "\"" <<endreq;
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
  primaryGenerator->init(primary, m_ppsvc);

  McParticleManager::getPointer()->addMcParticle(0,primary);   
  // Run geant4
  m_runManager->BeamOn(); 
    
  // If we set the pruneCal mode than we remove all the not TKR interacting
  // particles that does not touch the TKR zone (z>0)
  if (m_mcTreeMode=="pruneCal")
    McParticleManager::getPointer()->pruneCal();

  // Save the McParticle hierarchy in the TDS
  McParticleManager::getPointer()->save();
    
  // set up display of trajectories
  DisplayManager* dm = DisplayManager::instance();
  if(dm !=0) {   
    for( int i = 0; i< m_runManager->getNumberOfTrajectories(); ++i){
      std::auto_ptr<std::vector<Hep3Vector> > points = 
        m_runManager->getTrajectoryPoints(i);
      dm->addTrack(*(points.get()), m_runManager->getTrajectoryCharge(i), i==0);
    }
  }

  // we save the trajectories if the m_saveTrajectories is true
  if (m_saveTrajectories)
    {
      Event::McTrajectoryList* traj = new Event::McTrajectoryList();  
      eventSvc()->registerObject("Event/MC/TrajectoryCol",traj);

      for( int j = 0; j< m_runManager->getNumberOfTrajectories(); ++j){
        std::auto_ptr<std::vector<Hep3Vector> > points = 
          m_runManager->getTrajectoryPoints(j);

        Event::McParticle* part = 0;
      
        if (McParticleManager::getPointer()->size() > m_runManager->getTrajectoryTrackId(j))
          {
            part = 
              McParticleManager::getPointer()->getMcParticle(m_runManager->getTrajectoryTrackId(j));
          }
        
        if(part)
          {
            Event::McTrajectory* tr = new Event::McTrajectory();
            tr->addPoints(*(points.get()));
            tr->setMcParticle(part);
            
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





