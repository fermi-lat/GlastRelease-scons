// File and Version Information:
// $Header$
//
// Description: 
// This class manages the Geant4 main loop and calls; since we don't need event
// loop management, graphics interfaces or data persistency with standard Geant4
// features (we are using Gaudi for that, this is a stripped version of a true
// G4RunManager class. Most of its methods and members are hinerited from the
// Geant4 one; since they are not supposed to be used directly by any GLAST
// client (since they are automatically used by Geant4 internal architecture)
// they are not documented. For their meaning please see the standard Geant4
// documentation
//
// Author(s):
//      R.Giannitrapani

#include "G4Timer.hh"

#include "DetectorConstruction.h"
#include "RunManager.h"
#include "UIsession.h"
#include "PhysicsList.h"
#include "PrimaryGeneratorAction.h"
#include "TrackingAction.h"

#include "Randomize.hh"
#include "G4Run.hh"
#include "G4RunMessenger.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"

#include "G4UserRunAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeometryManager.hh"
#include "G4SDManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4VPersistencyManager.hh"
#include "G4UImanager.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessTable.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryPoint.hh"

#include "G4ios.hh"
#include "g4std/strstream"

#include <memory>

RunManager* RunManager::fRunManager = NULL;

RunManager* RunManager::GetRunManager()
{ return fRunManager; }

RunManager::RunManager(IGlastDetSvc* gds, IDataProviderSvc* esv,
                       std::string geometryMode)
  :userDetector(NULL),physicsList(NULL),
   userPrimaryGeneratorAction(NULL),
   currentRun(NULL),currentEvent(NULL),
   geometryInitialized(false),physicsInitialized(false),
   cutoffInitialized(false),
   geometryNeedsToBeClosed(true),initializedAtLeastOnce(false),
   runAborted(false),
   geometryToBeOptimized(true),verboseLevel(0),DCtable(NULL),runIDCounter(0),
   storeRandomNumberStatus(0)
{
  if(fRunManager)
    { G4Exception("RunManager constructed twice."); }

  fRunManager = this;

  // This dummy session is needed later to silent G4
  session = new UIsession;

  // The event manager of G4
  eventManager = new G4EventManager();
  // The timer of G4
  timer = new G4Timer();

  // Set the TrackingAction to track the McParticle
  eventManager->SetUserAction(new TrackingAction);

  // Various G4 messenger needed
  G4ParticleTable::GetParticleTable()->CreateMessenger();
  G4ProcessTable::GetProcessTable()->CreateMessenger();
  randomNumberStatusDir = "./";

  // The user stuff
  userDetector = new DetectorConstruction(gds,esv, geometryMode);
  physicsList = new PhysicsList;
  userPrimaryGeneratorAction = new PrimaryGeneratorAction;
}

RunManager::~RunManager()
{
  G4StateManager* pStateManager = G4StateManager::GetStateManager();
  pStateManager->SetNewState(Quit);

  delete session;

  delete timer;

  physicsList->RemoveProcessManager();
  G4ParticleTable::GetParticleTable()->DeleteMessenger();
  G4ProcessTable::GetProcessTable()->DeleteMessenger();

  if(userDetector)
    {
      delete userDetector;
    }
  if(physicsList)
    {
      delete physicsList;
    }
  if(userPrimaryGeneratorAction)
    {
      delete userPrimaryGeneratorAction;
    }
  G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
  if(fSDM)
    {
      delete fSDM;
    }
  delete eventManager;

  G4UImanager* pUImanager = G4UImanager::GetUIpointer();
  {
    if(pUImanager) delete pUImanager;
  }
  if(pStateManager) 
    {
      delete pStateManager;
    }
}

void RunManager::BeamOn()
{
  G4bool cond = ConfirmBeamOnCondition();
  if(cond)
    {    
      G4StateManager* stateManager = G4StateManager::GetStateManager();
      RunInitialization();
      stateManager->SetNewState(EventProc);
      currentEvent = GenerateEvent(0);
      eventManager->ProcessOneEvent(currentEvent);
      stateManager->SetNewState(GeomClosed);
      RunTermination();
    }
}

G4Event* RunManager::getCurrentEvent()const
{
  return currentEvent;
}

G4bool RunManager::ConfirmBeamOnCondition()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();

  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=PreInit && currentState!=Idle)
    {
      G4cerr << "Illegal application state - BeamOn() ignored." << G4endl;
      return false;
    }

  if(!initializedAtLeastOnce)
    {
      G4cerr << " Geant4 kernel should be initialized" << G4endl;
      G4cerr << "before the first BeamOn(). - BeamOn ignored." << G4endl;
      return false;
    }

  if(!(geometryInitialized && physicsInitialized && cutoffInitialized)) 
    {
      if(verboseLevel>0)
        {
          G4cout << "Start re-initialization because " << G4endl;
          if(!geometryInitialized) G4cout << "  Geometry" << G4endl;
          if(!physicsInitialized)  G4cout << "  Physics processes" << G4endl;
          if(!cutoffInitialized)   G4cout << "  Cutoff" << G4endl;
          G4cout << "has been modified since last Run." << G4endl;
        }
      Initialize();
    }

  return true;
}

void RunManager::RunInitialization()
{
  if (currentEvent) delete currentEvent;
  currentRun = new G4Run();
  currentRun->SetRunID(runIDCounter);

  currentRun->SetDCtable(DCtable);
  G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
  if(fSDM)
    { currentRun->SetHCtable(fSDM->GetHCtable()); }
  
  if(geometryNeedsToBeClosed)
    {
      if(verboseLevel>1) G4cout << "Start closing geometry." << G4endl;
      G4GeometryManager* geomManager = G4GeometryManager::GetInstance();
      geomManager->OpenGeometry();
      geomManager->CloseGeometry(geometryToBeOptimized);
      geometryNeedsToBeClosed = false;
    }
  G4StateManager* stateManager = G4StateManager::GetStateManager();
  stateManager->SetNewState(GeomClosed);

  runAborted = false;

  if(storeRandomNumberStatus==1 || storeRandomNumberStatus==-1)
    StoreRandomNumberStatus();
  
  if(verboseLevel>0) G4cout << "Start Run processing." << G4endl;
}

G4Event* RunManager::GenerateEvent(G4int i_event)
{
  if(!userPrimaryGeneratorAction)
    {
      G4Exception
        ("RunManager::BeamOn - G4VUserPrimaryGeneratorAction is not defined.");
    }

  G4Event* anEvent = new G4Event(i_event);

  if(storeRandomNumberStatus==2 || storeRandomNumberStatus==-2) 
    StoreRandomNumberStatus(anEvent->GetEventID());

  userPrimaryGeneratorAction->GeneratePrimaries(anEvent);
  return anEvent;
}

void RunManager::RunTermination()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();

  delete currentRun;
  currentRun = NULL;

  stateManager->SetNewState(Idle);
}

void RunManager::Initialize()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=PreInit && currentState!=Idle)
    {
      G4cerr << "Illegal application state - "
             << "RunManager::Initialize() ignored." << G4endl;
      return;
    }

  /// Set a dummy session to silent G4 intitialization messages on screen
  G4UImanager* pUImanager = G4UImanager::GetUIpointer();
  pUImanager->SetCoutDestination(session);

  stateManager->SetNewState(Init);
  if(!geometryInitialized) InitializeGeometry();
  if(!physicsInitialized) InitializePhysics();
  if(!cutoffInitialized) InitializeCutOff();
  stateManager->SetNewState(Idle);
  if(!initializedAtLeastOnce) initializedAtLeastOnce = true;

  pUImanager->SetCoutDestination(new G4UIsession);
}

void RunManager::InitializeGeometry()
{
  if(!userDetector)
    {
      G4Exception
        ("G4VUserDetectorConstruction is not defined.");
    }

  if(verboseLevel>1) G4cout << "userDetector->Construct() start." << G4endl;
  DefineWorldVolume(userDetector->Construct());
  geometryInitialized = true;
}

void RunManager::InitializePhysics()
{
  if(physicsList)
    {
      if(verboseLevel>1) G4cout << "physicsList->Construct() start." << G4endl;
      physicsList->SetVerboseLevel(0);
      physicsList->Construct();
    }
  else
    {
      G4Exception("G4VUserPhysicsList is not defined");
    }
  physicsInitialized = true;
}

void RunManager::InitializeCutOff()
{
  if(physicsList)
    {
      if(verboseLevel>1) G4cout << "physicsList->setCut() start." << G4endl;
      physicsList->SetCuts();
    }
  cutoffInitialized = true;
}
  
void RunManager::AbortRun()
{
  // This method is valid only for GeomClosed or EventProc state
  G4ApplicationState currentState = 
    G4StateManager::GetStateManager()->GetCurrentState();
  if(currentState==GeomClosed || currentState==EventProc)
    {
      runAborted = true;
      if(currentState==EventProc) eventManager->AbortCurrentEvent();
    }
  else
    {
      G4cerr << "Run is not in progress. AbortRun() ignored." << G4endl;
    }
}

void RunManager::DefineWorldVolume(G4VPhysicalVolume* worldVol)
{
  // set the world volume to the Navigator
  G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking()
    ->SetWorldVolume(worldVol);

  // Let VisManager know it
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager) pVVisManager->GeometryHasChanged();

  geometryNeedsToBeClosed = true;
}

void RunManager::StoreRandomNumberStatus(G4int eventID)
{
  G4String fileN = randomNumberStatusDir+"RandEngine";
  if(storeRandomNumberStatus>0 && currentRun != NULL)
    {
      char st[20];
      G4std::ostrstream os(st,20);
      os << currentRun->GetRunID() << '\0';
      fileN += "R";
      fileN += st;
    }
  if(storeRandomNumberStatus==2 && eventID>=0)
    {
      char st[20];
      G4std::ostrstream os(st,20);
      os << eventID << '\0';
      fileN += "E";
      fileN += st;
    }
  fileN += ".stat";
  HepRandom::saveEngineStatus(fileN);
}
  
void RunManager::RestoreRandomNumberStatus(G4String fileN)
{
  G4String fileNameWithDirectory;
  if(fileN.index("/")==G4std::string::npos)
    { fileNameWithDirectory = randomNumberStatusDir+fileN; }
  else
    { fileNameWithDirectory = fileN; }
  HepRandom::restoreEngineStatus(fileNameWithDirectory);
}


unsigned int RunManager::getNumberOfTrajectories()
{
  G4Event* currentEvent = getCurrentEvent();
  if (currentEvent->GetTrajectoryContainer())
    return(*((currentEvent)->GetTrajectoryContainer())).entries();
  else return 0;
}

int RunManager::getTrajectoryCharge(unsigned int i)
{
  if (i > getNumberOfTrajectories())
    return -99;
  
  G4Event* event = getCurrentEvent();
  if (event->GetTrajectoryContainer())
    {
      G4Trajectory* trajectory = 
        static_cast<G4Trajectory*>((*((event)->GetTrajectoryContainer()))[i]);
      
      return trajectory->GetCharge();
    }
  else return -99;
}

std::auto_ptr<std::vector<Hep3Vector> > RunManager::getTrajectoryPoints(unsigned int i)
{
  if (i > getNumberOfTrajectories())
    return std::auto_ptr<std::vector<Hep3Vector> >(0);
  
  G4Event* event = getCurrentEvent();
  if (event->GetTrajectoryContainer())
    {
      G4Trajectory* trajectory =
        static_cast<G4Trajectory*>((*((event)->GetTrajectoryContainer()))[i]);
      
      std::vector<Hep3Vector>* points = new std::vector<Hep3Vector>;
      for(unsigned int j=0;j<trajectory->GetPointEntries();j++)
        {  
          G4TrajectoryPoint* currentPoint = 
            static_cast<G4TrajectoryPoint*>(trajectory->GetPoint(j));
          points->push_back(currentPoint->GetPosition());
        }
      return std::auto_ptr<std::vector<Hep3Vector> >(points);
    }
  else return std::auto_ptr<std::vector<Hep3Vector> >(0);
}







