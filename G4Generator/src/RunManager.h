#ifndef RunManager_h
#define RunManager_h 1

// userAction classes
class UIsession;

class G4VUserDetectorConstruction;
class G4VUserPhysicsList;
class G4VModularPhysicsList;
class G4VUserPrimaryGeneratorAction;

class G4VPhysicalVolume;
class G4Timer;
class G4RunMessenger;
class G4DCtable;
class G4Run;

class Hep3Vector;

#include "G4Event.hh"

#include "G4EventManager.hh"
#include "globals.hh"
#include "g4std/vector"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include <memory>

/** 
 * @class RunManager
 *
 * @brief The Geant4 singleton manager
 *
 * This is a very draft version of a customized G4RunManager for Glast; it
 * overrides the normal eventloop of G4, letting Gaudi to manage it. Most of its
 * functionalities are duplicated from the Geant4 class and are not relevant for
 * the rest of the code. Some methods will be removed in the near future since
 * they are not used by GLAST. For these motivations this class is not
 * completely documented. Please see the Geant4 toolkit for more information.
 *
 * Please note that overriding the RunManager of Geant4, while has been quite
 * easy and fast to do and is providing us with the expected result, is risky
 * since we are not garanteeded against future releases of Geant4. 
 *
 * @author R.Giannitrapani
 *    
 * $Header$
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
     The constructor needs a pointer to the abstract interface of the
     GlastDetSvc and to the DataProviderSvc. It gets also the mode for the
     geometry level of details
  */
  RunManager(IGlastDetSvc* gds, IDataProviderSvc* esv, 
             std::string geometryMode);
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

  /// This method return a pointer to the current G4Event
  G4Event* getCurrentEvent() const;
  
  /// This method return the number of trajectories stored in the currentEvent
  unsigned int getNumberOfTrajectories();

  /// This method return the vector of points of a trajectory
  int getTrajectoryCharge(unsigned int);

  /// This method return the vector of points of a trajectory
  std::auto_ptr<std::vector<Hep3Vector> > getTrajectoryPoints(unsigned int);

  
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

  //  G4VUserPhysicsList * physicsList;

  G4VModularPhysicsList * physicsList;

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
  inline const G4VModularPhysicsList* GetUserPhysicsList() const
    { return physicsList; }
  inline const G4VUserPrimaryGeneratorAction* GetUserPrimaryGeneratorAction() 
    const { return userPrimaryGeneratorAction; }
  
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
  /// Returns the pointer to the current run. This method is available for
  /// Geant4 states of GeomClosed and EventProc.
  inline const G4Event* GetCurrentEvent() const
    { return currentEvent; }
  /// Returns the pointer to the current event. This method is available for
  /// EventProc state.
  inline void SetRunIDCounter(G4int i)
    { runIDCounter = i; }
  /// Set the run number counter. Initially, the counter is initialized to zero
  /// and incremented by one for every BeamOn().

  public:
    inline void SetDCtable(G4DCtable* DCtbl)
    { DCtable = DCtbl; }
};

#endif




