// $Header$
#include "G4Timer.hh"

#include "RunManager.h"
#include "UIsession.h"
#include "DetectorConstruction.h"
#include "PhysicsList.h"
#include "PrimaryGeneratorAction.h"

#include "Randomize.hh"
#include "G4Run.hh"
#include "G4RunMessenger.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
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

#include "G4ios.hh"
#include "g4std/strstream"


RunManager* RunManager::fRunManager = NULL;

RunManager* RunManager::GetRunManager()
{ return fRunManager; }

RunManager::RunManager(std::string topvol, std::string visitorMode)
  :userDetector(NULL),physicsList(NULL),
   userPrimaryGeneratorAction(NULL),
   currentRun(NULL),currentEvent(NULL),
   geometryInitialized(false),physicsInitialized(false),cutoffInitialized(false),
   geometryNeedsToBeClosed(true),initializedAtLeastOnce(false),
   runAborted(false),
   geometryToBeOptimized(true),verboseLevel(0),DCtable(NULL),runIDCounter(0),
   storeRandomNumberStatus(0),
   m_topvol(topvol),
   m_visitorMode(visitorMode)
{
  if(fRunManager)
  { G4Exception("RunManager constructed twice."); }
  //G4UnitDefinition::BuildUnitsTable();
  fRunManager = this;

  session = new UIsession;

  eventManager = new G4EventManager();
  timer = new G4Timer();

  G4ParticleTable::GetParticleTable()->CreateMessenger();
  G4ProcessTable::GetProcessTable()->CreateMessenger();
  randomNumberStatusDir = "./";

  // The user stuff
  userDetector = new DetectorConstruction(m_topvol, m_visitorMode);
  physicsList = new PhysicsList;
  userPrimaryGeneratorAction = new PrimaryGeneratorAction;

  G4UImanager* pUImanager = G4UImanager::GetUIpointer();
  pUImanager->SetCoutDestination(session);
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
    currentEvent = NULL;

    RunTermination();
  }
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

  if(storeRandomNumberStatus==1 || storeRandomNumberStatus==-1) StoreRandomNumberStatus();
  
  if(verboseLevel>0) G4cout << "Start Run processing." << G4endl;
}

void RunManager::DoEventLoop(G4int n_event,const char* macroFile,G4int n_select)
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();

  if(verboseLevel>0) 
  { timer->Start(); }

  G4String msg;
  if(macroFile!=NULL)
  { 
    if(n_select<0) n_select = n_event;
    msg = "/control/execute ";
    msg += macroFile;
  }
  else
  { n_select = -1; }

  G4int i_event;
  for( i_event=0; i_event<n_event; i_event++ )
  {

  }

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

void RunManager::AnalyzeEvent(G4Event* anEvent)
{
}

void RunManager::RunTermination()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();

  delete currentRun;
  currentRun = NULL;

  stateManager->SetNewState(Idle);
}

void RunManager::StackPreviousEvent(G4Event* anEvent)
{
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

  stateManager->SetNewState(Init);
  if(!geometryInitialized) InitializeGeometry();
  if(!physicsInitialized) InitializePhysics();
  if(!cutoffInitialized) InitializeCutOff();
  stateManager->SetNewState(Idle);
  if(!initializedAtLeastOnce) initializedAtLeastOnce = true;
}

void RunManager::InitializeGeometry()
{
  if(!userDetector)
  {
    G4Exception
    ("RunManager::InitializeGeometry - G4VUserDetectorConstruction is not defined.");
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








