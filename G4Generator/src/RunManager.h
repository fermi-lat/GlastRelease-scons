// $Header$

/** Very draft version of a customized G4RunManager for Glast;
    It will override the normal eventloop of G4, letting Gaudi
    to manage it
*/

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

class RunManager
{
 public: // with description
  static RunManager* GetRunManager();
  //  Static method which returns the singleton pointer of RunManager or
  // its derived class.
  
 private:
  static RunManager* fRunManager;
  
 public: // with description
  RunManager();
  virtual ~RunManager();
  //  The constructor and the destructor. The user must construct this class
  // object at the beginning of his/her main() and must delete it at the 
  // bottom of the main().

  public: // with description
    virtual void BeamOn();

    virtual void Initialize();
    //  This method invokes all the necessary initialization procedures for an event
    // loop. This method must be invoked at the Geant4 state of PreInit or Idle. The
    // state will be changed to Init during the initialization procedures and then
    // changed to Idle.
    //  This method invokes three protected methods, InitializeGeometry(), 
    // InitializePhysics(), and InitializeCutOff().
    //  After some event loops, the user can invoke this method once again. It is
    // required if the user changes geometry, physics process, and/or cut off value.
    // If the user forget the second invokation, RunManager will invoke BeamOn()
    // method will invoke this method. (Note that this feature is not valid for the
    // first initialization.)
    virtual void DefineWorldVolume(G4VPhysicalVolume * worldVol);
    //  Usually, this method is invoked from InitializeGeometry() protected method
    // of this class. But, in case all of geometry has already created and kept in
    // the ODBMS, the pointer to the world physical volume can be set by this method.
    virtual void AbortRun();
    //  This method safely aborts the current event loop even if an event is in progress.
    // This method is available for Geant4 states of GeomClosed and EventProc. The state
    // will be changed to Idle, so that another event loop can be done.

  protected: // with description

    virtual void InitializeGeometry();
    virtual void InitializePhysics();
    virtual void InitializeCutOff();
    //  These three protected methods are invoked from Initialize() method for the 
    // initializations of geometry, physics processes, and cut off. The user's concrete
    // G4VUserDetectorConstruction class will be accessed from InitializeGeometry() and
    // G4VUserPhysicsList class will be accessed from other two methods.

    virtual G4bool ConfirmBeamOnCondition();
    virtual void RunInitialization();
    virtual void DoEventLoop(G4int n_event,const char* macroFile=NULL,G4int n_select=-1);
    virtual void RunTermination();
    //  These four protected methods are invoked from BeamOn() method. These four methods
    // are invoked in this order.
    //  ConfirmBeamOnCondition() method checks if all the necessary initializations have
    // already done. If the condition is not satisfied, false is returned and the follwing
    // three methods will be skipped.
    //  RunInitialization() method initializes a run. For example, a G4Run class object 
    // is constructed in this method.
    //  DoEventLoop() method control an event loop. Arguments are same as BeamOn() method.
    // Inide the event loop, two following protected methods are invoked at the begining
    // and the end of each event. 
    //  RunTermination() method terminates a run processing. For example, a G4Run class
    // object is deleted in this class. If the user uses ODBMS and wants to store the
    // G4Run class object, he/she must override this method.

    virtual G4Event* GenerateEvent(G4int i_event);
    virtual void AnalyzeEvent(G4Event* anEvent);
    //  These two protected methods are invoked from DoEventLoop() method at the begining
    // and the end of each event processing.
    //  GenerateEvent() method constructs a G4Event class object and invoke the user's
    // G4VUserPrimaryGeneratorAction concrete class. If the user is using ODBMS and event 
    // objects have been created and stored in the data base, he/she must override this
    // method.
    //  AnalyzeEvent() stores an event to a data base if a concrete G4VPersistentManager
    // class is defined.

  protected:
    void StackPreviousEvent(G4Event* anEvent);

  protected:
    G4EventManager * eventManager;
    G4VUserDetectorConstruction * userDetector;
    G4VUserPhysicsList * physicsList;
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

    UIsession* session;
  public:
    virtual void StoreRandomNumberStatus(G4int eventID=-1);
    virtual void RestoreRandomNumberStatus(G4String fileN);

    //  These methods store respective user initialization and action classes.
    inline const G4VUserDetectorConstruction* GetUserDetectorConstruction() const
    { return userDetector; }
    inline const G4VUserPhysicsList* GetUserPhysicsList() const
    { return physicsList; }
    inline const G4VUserPrimaryGeneratorAction* GetUserPrimaryGeneratorAction() const
    { return userPrimaryGeneratorAction; }

  public:
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

  public: // with description
    inline void GeometryHasBeenModified()
    { geometryNeedsToBeClosed = true; }
    inline void CutOffHasBeenModified()
    { cutoffInitialized = false; }
    //  These two methods must be invoked (or equivalent UI commands can be used)
    // in case the user changes his/her detector geometry or cut off value(s) after
    // Initialize() metho has been invoked. Then, at the begining of the next BeamOn(),
    // all necessary re-initialization will be done.

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

  public: // with description
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

