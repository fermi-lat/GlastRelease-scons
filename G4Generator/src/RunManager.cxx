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

#include "RunManager.h"
#include "G4RunManagerKernel.hh"

#include "UIsession.h"
#include "PhysicsList.h"
#include "PrimaryGeneratorAction.h"
#include "TrackingAction.h"
#include "SteppingAction.h"
#include "G4Generator/IG4GeometrySvc.h"
#include "Randomize.hh"

#include "G4StateManager.hh"
#include "G4ApplicationState.hh"
#include "Randomize.hh"
#include "G4Run.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4TransportationManager.hh"
#include "G4VPersistencyManager.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessTable.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryPoint.hh"

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4ios.hh"
#include <strstream>

#include "src/G4GenExceptionHandler.h"

#include "G4PhysListFactory.hh"

RunManager* RunManager::fRunManager = 0;

RunManager* RunManager::GetRunManager()
{ return fRunManager; }

RunManager::RunManager(std::ostream& log, 
                       double defaultCutValue, 
                       double defaultTkrCutValue,
                       double defaultCalCutValue,
                       std::string& physics_choice, 
                       std::string& physics_table,
                       std::string&  physics_dir,
                       GlastMS::MultipleScatteringFactory& msfactory,
                       GlastMS::EnergyLossFactory& eLossFactory,
					   IG4GeometrySvc* gsv)
:m_log(log),
 physicsList(NULL),
 userPrimaryGeneratorAction(NULL),
 geometryInitialized(false),
 physicsInitialized(false),
 geometryNeedsToBeClosed(true),
 runAborted(false),
 initializedAtLeastOnce(false),
 geometryToBeOptimized(true),
 currentEvent(0),
 runIDCounter(0),
 verboseLevel(0),
 defaultCut(defaultCutValue),
 TkrCutValue(defaultTkrCutValue),
 CalCutValue(defaultCalCutValue),
 storeRandomNumberStatus(0),
 m_gsv(gsv)
{
  // Create a (singleton) Run Manager
  if(fRunManager)
  { G4Exception("RunManager constructed twice."); }
  fRunManager = this;

  // This dummy session is needed later to silent G4
  session = new UIsession;

  // Now create a G4RunManagerKernel
  kernel       = new G4RunManagerKernel();
  eventManager = kernel->GetEventManager();

  // G4 timer
  timer = new G4Timer();

  // Set the TrackingAction to track the McParticle
  //eventManager->SetUserAction(new TrackingAction(gsv));
  //eventManager->SetUserAction(new SteppingAction());

  // G4 Messenger stuff
  G4ParticleTable::GetParticleTable()->CreateMessenger();
  G4ProcessTable::GetProcessTable()->CreateMessenger();
  randomNumberStatusDir = "./";


  


  // The user stuff

  //  physicsList = new PhysicsList(defaultCutValue, physics_choice, physics_table, physics_dir, msfactory, eLossFactory);
  
  G4PhysListFactory* myFactory = new G4PhysListFactory();
  
  physicsList = myFactory -> GetReferencePhysList(G4String(physics_choice));

  //userPrimaryGeneratorAction = new PrimaryGeneratorAction;



}

RunManager::~RunManager()
{
  delete timer;
  G4ParticleTable::GetParticleTable()->DeleteMessenger();
  G4ProcessTable::GetProcessTable()->DeleteMessenger();
  if(physicsList)
  {
    delete physicsList;
    if(verboseLevel>1) G4cout << "UserPhysicsList deleted." << G4endl;
  }
  if(userPrimaryGeneratorAction)
  {
    delete userPrimaryGeneratorAction;
    if(verboseLevel>1) G4cout << "UserPrimaryGenerator deleted." << G4endl;
  }

  delete kernel;

  if(verboseLevel>1) G4cout << "RunManager is deleting." << G4endl;
}

void RunManager::BeamOn()
{
  G4bool cond = ConfirmBeamOnCondition();
  if(cond)
    {    
//      G4StateManager* stateManager = G4StateManager::GetStateManager();
      RunInitialization();
//      stateManager->SetNewState(G4State_EventProc);
      currentEvent = GenerateEvent(0);
      eventManager->ProcessOneEvent(currentEvent);
//      stateManager->SetNewState(G4State_GeomClosed);
      RunTermination();
    }
}

G4bool RunManager::ConfirmBeamOnCondition()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();

  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=G4State_PreInit && currentState!=G4State_Idle)
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
  if(!(kernel->RunInitialization())) return;
  currentRun = new G4Run();

  currentRun->SetRunID(runIDCounter);

  currentRun->SetDCtable(DCtable);
  G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
  if(fSDM)
  { currentRun->SetHCtable(fSDM->GetHCtable()); }

  // Insert here something to clean up the current event 
  // which appears to have been a memory leak in our code
  if (currentEvent)
  {
      delete currentEvent;
      currentEvent = 0;
  }

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
  delete currentRun;
  currentRun = 0;

  kernel->RunTermination();
}


void RunManager::Initialize()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=G4State_PreInit && currentState!=G4State_Idle)
  {
    G4cerr << "Illegal application state - "
         << "RunManager::Initialize() ignored." << G4endl;
    return;
  }

  /// Set a dummy session to silent G4 intitialization messages on screen
  G4UImanager* pUImanager = G4UImanager::GetUIpointer();
  
  pUImanager->SetCoutDestination(session);
  
  if(!geometryInitialized) InitializeGeometry();
  if(!physicsInitialized) InitializePhysics();
  if(!cutoffInitialized) InitializeCutOff();
  if(!initializedAtLeastOnce) initializedAtLeastOnce = true;

  userPrimaryGeneratorAction = new PrimaryGeneratorAction;

  // Set the TrackingAction to track the McParticle
  eventManager->SetUserAction(new TrackingAction(m_gsv));
  eventManager->SetUserAction(new SteppingAction(m_gsv));

  pUImanager->SetCoutDestination(new G4UIsession);

  // Set the exception handler
  m_ExceptionHandler = new G4GenExceptionHandler();
}

void RunManager::InitializeGeometry()
{
  // Check the world volume to the Navigator
  if (G4VPhysicalVolume* currentWorld = G4TransportationManager::GetTransportationManager()
              ->GetNavigatorForTracking()->GetWorldVolume())
  {
    kernel->DefineWorldVolume(currentWorld, false);

    // Check to see if the geometry has already been associated to a region
    if (!currentWorld->GetLogicalVolume()->GetRegion())
    {
        // If so, this is an error
        G4cerr << "The world volume has a user-defined region <" 
               << currentWorld->GetLogicalVolume()->GetRegion()->GetName()
               << ">." << G4endl;
    }
    else
    {
        // Define the default region for the world geometry
//        G4Region* defaultRegion = new G4Region("DefaultRegionForGLAST");
//        defaultRegion->SetProductionCuts(G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());
        
        // set the default region to the world
//        G4LogicalVolume* worldLog = currentWorld->GetLogicalVolume();
//        worldLog->SetRegion(defaultRegion);
//        defaultRegion->AddRootLogicalVolume(worldLog);

        // Now want to define special regions for the track and the cal. To do this
        // we need to search through the logical volumes to find the mother volumes for
        // the track and the cal. 
        G4LogicalVolumeStore* Store    = G4LogicalVolumeStore::GetInstance();
        int                   nVolumes = Store->size();
        G4LogicalVolumeStore::iterator volIter = Store->begin();

        bool needTkr = true;
        bool needCal = true;

        // Loop over volumes until end or we are finished
        while(volIter < Store->end() && (needTkr || needCal))
        {
            G4LogicalVolume* logVol  = *volIter++;

            // Is it the Tracker mother volume?
            if (logVol->GetName() == "oneTKR")
            {
                // Define a special region for the Tracker
                G4Region*         tkrRegion = new G4Region("TrackerRegion");
                G4ProductionCuts* tkrCuts   = new G4ProductionCuts();
                tkrCuts->SetProductionCut(TkrCutValue);
                tkrRegion->SetProductionCuts(tkrCuts);

                logVol->SetRegion(tkrRegion);
                tkrRegion->AddRootLogicalVolume(logVol);
                needTkr = false;
            }

            // Is it the Cal mother volume?
            if (logVol->GetName() == "oneCAL")
            {
                // Define a region for the Cal
                G4Region*         calRegion = new G4Region("CalorimeterRegion");
                G4ProductionCuts* calCuts   = new G4ProductionCuts();
                calCuts->SetProductionCut(CalCutValue);
                calRegion->SetProductionCuts(calCuts);

                logVol->SetRegion(calRegion);
                calRegion->AddRootLogicalVolume(logVol);
                needCal = false;
            }
        }
    }

	// Let VisManager know it
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	if(pVVisManager) pVVisManager->GeometryHasChanged();

	geometryNeedsToBeClosed = true;
  }
  else
  {
    G4cerr << "Cannot find a valid geometry! " << G4endl; 
  }
  geometryInitialized = true;
}

void RunManager::InitializePhysics()
{
  if(physicsList)
  {
    if(verboseLevel>1) G4cout << "physicsList->Construct() start." << G4endl;
    //kernel->InitializePhysics(physicsList);   // Use this for G4 v6
    kernel->SetPhysics(physicsList);        // Use this for G4 v7
    kernel->InitializePhysics();
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

void RunManager::AbortRun(G4bool softAbort)
{
  // This method is valid only for GeomClosed or EventProc state
  G4ApplicationState currentState = 
    G4StateManager::GetStateManager()->GetCurrentState();
  if(currentState==G4State_GeomClosed || currentState==G4State_EventProc)
  {
    runAborted = true;
    if(currentState==G4State_EventProc && !softAbort)
    {
      currentEvent->SetEventAborted();
      eventManager->AbortCurrentEvent();
    }
  }
  else
  {
    G4cerr << "Run is not in progress. AbortRun() ignored." << G4endl;
  }
}

void RunManager::AbortEvent()
{
  // This method is valid only for EventProc state
  G4ApplicationState currentState = 
    G4StateManager::GetStateManager()->GetCurrentState();
  if(currentState==G4State_EventProc)
  {
    currentEvent->SetEventAborted();
    eventManager->AbortCurrentEvent();
  }
  else
  {
    G4cerr << "Event is not in progress. AbortEevnt() ignored." << G4endl;
  }
}

void RunManager::StoreRandomNumberStatus(G4int eventID)
{
  G4String fileN = randomNumberStatusDir+"RandEngine";
  if(storeRandomNumberStatus>0 && currentRun != NULL)
    {
      char st[20];
      std::ostrstream os(st,20);
      os << currentRun->GetRunID() << '\0';
      fileN += "R";
      fileN += st;
    }
  if(storeRandomNumberStatus==2 && eventID>=0)
    {
      char st[20];
      std::ostrstream os(st,20);
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
  if(fileN.index("/")==std::string::npos)
  { fileNameWithDirectory = randomNumberStatusDir+fileN; }
  else
  { fileNameWithDirectory = fileN; }
  
  HepRandom::restoreEngineStatus(fileNameWithDirectory);
  if(verboseLevel>0) G4cout << "RandomNumberEngineStatus restored from file: "
         << fileNameWithDirectory << G4endl;
  HepRandom::showEngineStatus();	 
}

void RunManager::DumpRegion(G4String rname) const
{
  kernel->UpdateRegion();
  kernel->DumpRegion(rname);
}

void RunManager::DumpRegion(G4Region* region) const
{
  kernel->UpdateRegion();
  kernel->DumpRegion(region);
}


unsigned int RunManager::getNumberOfTrajectories()
{
  const G4Event* currentEvent = GetCurrentEvent();
  if (currentEvent->GetTrajectoryContainer())
    return(*((currentEvent)->GetTrajectoryContainer())).entries();
  else return 0;
}

int RunManager::getTrajectoryCharge(unsigned int i)
{
  if (i > getNumberOfTrajectories())
    return -99;
  
  const G4Event* event = GetCurrentEvent();
  if (event->GetTrajectoryContainer())
    {
      G4Trajectory* trajectory = 
        static_cast<G4Trajectory*>((*((event)->GetTrajectoryContainer()))[i]);
      
      return (int)trajectory->GetCharge();
    }
  else return -99;
}

std::auto_ptr<std::vector<Hep3Vector> > RunManager::getTrajectoryPoints(unsigned int i)
{
  if (i > getNumberOfTrajectories())
    return std::auto_ptr<std::vector<Hep3Vector> >(0);
  
  const G4Event* event = GetCurrentEvent();
  if (event->GetTrajectoryContainer())
    {
      G4Trajectory* trajectory =
        static_cast<G4Trajectory*>((*((event)->GetTrajectoryContainer()))[i]);
      
      std::vector<Hep3Vector>* points = new std::vector<Hep3Vector>;
      for(int j=0;j<trajectory->GetPointEntries();j++)
        {  
          G4TrajectoryPoint* currentPoint = 
            static_cast<G4TrajectoryPoint*>(trajectory->GetPoint(j));
          points->push_back(currentPoint->GetPosition());
        }
      return std::auto_ptr<std::vector<Hep3Vector> >(points);
    }
  else return std::auto_ptr<std::vector<Hep3Vector> >(0);
}

int RunManager::getTrajectoryTrackId(unsigned int i)
{
  if (i > getNumberOfTrajectories())
    return -1;
  
  const G4Event* event = GetCurrentEvent();
  if (event->GetTrajectoryContainer())
    {
      G4Trajectory* trajectory =
        static_cast<G4Trajectory*>((*((event)->GetTrajectoryContainer()))[i]);
      
      return trajectory->GetTrackID();
    }
  else return -1;
}


