// $Header$

#ifndef RunManager_h
#define RunManager_h 1

// userAction classes
class UIsession;

class G4VUserDetectorConstruction;
class G4VUserPhysicsList;
class G4VUserPrimaryGeneratorAction;

class G4VPhysicalVolume;
class G4Timer;
class G4RunMessenger;
class G4DCtable;
class G4Run;
#include "G4Event.hh"

#include "G4EventManager.hh"
#include "globals.hh"
#include "g4std/vector"

/** Very draft version of a customized G4RunManager for Glast;
    It will override the normal eventloop of G4, letting Gaudi
    to manage it.
*/
class RunManager
{
 public: 
  ///  Static method which returns the singleton pointer of RunManager
  static RunManager* GetRunManager();
  
 private:
  /// This is the pointer to the RunManager
  static RunManager* fRunManager;
  /// This is the top volume used for geometry
  std::string m_topvol;

  /// mode to apply to start the visitor, like "fastmc"
  std::string m_visitorMode;
  
 public: 
  /** 
     The constructor, with the top volume to be visited for geometry and the 
     mode to be used in the visit
  */
  RunManager(std::string topvol, std::string visitorMode);
  virtual ~RunManager();

 public: 
  /// This is the method to be invoked to start the simulation
  virtual void BeamOn();

  /// Set up the G4 stuff
  virtual void Initialize();

  /// Not used
  virtual void DefineWorldVolume(G4VPhysicalVolume * worldVol);

  /// Not used
  virtual void AbortRun();

  G4Event* getCurrentEvent() const;
  
 protected: 

  ///  These three protected methods are invoked from Initialize() method
  virtual void InitializeGeometry();
  virtual void InitializePhysics();
  virtual void InitializeCutOff();

  /// Confirm that G4 is in the right state to start an Event
  virtual G4bool ConfirmBeamOnCondition();
  /// Init and termination of a RUN .. to be eliminated?
  virtual void RunInitialization();
  virtual void RunTermination();

  /// Method to generate the Event .. the argument must be eliminated
  virtual G4Event* GenerateEvent(G4int i_event);

  protected:
  /// The Event manager of G4
  G4EventManager * eventManager;
  /// The detector construction class
  G4VUserDetectorConstruction * userDetector;
  /// The physics list class
  G4VUserPhysicsList * physicsList;
  /// The primary generator class
  G4VUserPrimaryGeneratorAction * userPrimaryGeneratorAction;

 protected:
  G4bool geometryInitialized;
  G4bool physicsInitialized;
  G4bool cutoffInitialized;
  G4bool geometryNeedsToBeClosed;
  G4bool runAborted;
  G4bool initializedAtLeastOnce;
  G4bool geometryToBeOptimized;

  G4int runIDCounter;
  G4int verboseLevel;
  G4Timer * timer;
  G4DCtable* DCtable;

  G4Run* currentRun;
  G4Event* currentEvent;

  G4int storeRandomNumberStatus;
  G4String randomNumberStatusDir;

  /// This is a dummy session to be used to silent G4
  UIsession* session;

 public:
  virtual void StoreRandomNumberStatus(G4int eventID=-1);
  virtual void RestoreRandomNumberStatus(G4String fileN);

  /// Return methods for the 3 user classes
  inline const G4VUserDetectorConstruction* GetUserDetectorConstruction() const
    { return userDetector; }
  inline const G4VUserPhysicsList* GetUserPhysicsList() const
    { return physicsList; }
  inline const G4VUserPrimaryGeneratorAction* GetUserPrimaryGeneratorAction() const
    { return userPrimaryGeneratorAction; }
  
  public:

  /// Stuff for random distributions ... ???
  inline void SetRandomNumberStore(G4int i)
    { storeRandomNumberStatus = i; }
  inline G4int GetRandomNumberStore() const
    { return storeRandomNumberStatus; }
  inline void SetRandomNumberStoreDir(G4String dir)
    { 
      G4String dirStr = dir;
      if( dirStr(dirStr.length()-1) != '/' ) dirStr += "/";
      randomNumberStatusDir = dirStr;
    }
  inline G4String GetRandomNumberStoreDir() const
    { return randomNumberStatusDir; }
  
 public: 
  /// This methods are necessary ???
  inline void GeometryHasBeenModified()
    { geometryNeedsToBeClosed = true; }
  inline void CutOffHasBeenModified()
    { cutoffInitialized = false; }
  

  public:
  inline void SetVerboseLevel(G4int vl)
    { verboseLevel = vl; }
  inline G4int GetVerboseLevel() const
    { return verboseLevel; }

  inline void SetGeometryToBeOptimized(G4bool vl)
    { 
      if(geometryToBeOptimized != vl)
	{
	  geometryToBeOptimized = vl;
	  geometryNeedsToBeClosed = true;
	}
    }
  inline G4bool GetGeometryToBeOptimized()
    { return geometryToBeOptimized; }
  
 public: 
  inline const G4Run* GetCurrentRun() const
    { return currentRun; }
  //  Returns the pointer to the current run. This method is available for Geant4
  // states of GeomClosed and EventProc.
  inline const G4Event* GetCurrentEvent() const
    { return currentEvent; }
  //  Returns the pointer to the current event. This method is available for EventProc
  // state.
  inline void SetRunIDCounter(G4int i)
    { runIDCounter = i; }
  //  Set the run number counter. Initially, the counter is initialized to zero and
  // incremented by one for every BeamOn().

  public:
    inline void SetDCtable(G4DCtable* DCtbl)
    { DCtable = DCtbl; }
};

#endif




