
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
//
//
// 

// class description:
//
//      This is a class for run control in GEANT4
// 
//     User must provide his own classes derived from the following
//     three abstract classes and register them to the RunManager. 
//        G4VUserPhysicsList                - Particle types and Processes 
//        G4VUserPrimaryGeneratorAction     - Event Generator selection
// 
//     In addition to the above mandatory classes, user can easily 
//     customize of the default functionality of GEANT4 simulation
//     by making his own classes derived from the following 5 user
//     action classes. 
//         G4UserEventAction                 - Actions for each Event
//         G4UserStackingAction              - Tracks Stacking selection
//         G4UserTrackingAction              - Actions for each Track
//         G4UserSteppingAction              - Actions for each Step
//     
//     RunManager is the only manager class in Geant4 kernel which 
//     the user MUST construct an object by him/herself in the main(). 
//     Also, RunManager is the only manager class in Geant4 kernel
//     which the user CAN derive it to costomize the behavior of the
//     run control. For this case, user should use protected methods
//     provided in this class for procedures he/she does not want to
//     change.
//
//     RunManager or the derived class of it MUST be a singleton.
//     The user MUST NOT construct more than one object even if there
//     are two different concrete implementations.
//
//     RunManager controls all of state changes. See G4ApplicationState.hh
//     in intercoms category for the meanings of each state.
//

#ifndef RunManager_h
#define RunManager_h 1

// userAction classes
class UIsession;

class G4VUserPhysicsList;
class G4VModularPhysicsList;
class G4VUserPrimaryGeneratorAction;

class G4VPhysicalVolume;
class G4Region;
class G4Timer;
class G4DCtable;
class G4Run;

namespace CLHEP {class Hep3Vector;}
using namespace CLHEP;

class IG4GeometrySvc;

class G4GenExceptionHandler;

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "G4RunManagerKernel.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "globals.hh"

#include <vector>
#include <memory>

namespace GlastMS { class MultipleScatteringFactory; class EnergyLossFactory;}

class RunManager
{
  public: // with description
    ///  Static method which returns the singleton pointer of RunManager
    static RunManager* GetRunManager();

  private:
    /// This is the pointer to the RunManager
    static RunManager* fRunManager;

    /// log for summary output
    std::ostream& m_log;

    /// Exception handler for G4 code
    G4GenExceptionHandler* m_ExceptionHandler;

  public: 
    /** 
       The constructor needs a pointer to the abstract interface of the
       GlastDetSvc and to the DataProviderSvc. It gets also the mode for the
       geometry level of details
    */
    RunManager(std::ostream& log, 
             double defaultCutValue, 
             double defaultTkrCutValue,
             double defaultCalCutValue,
             std::string& physics_choice, 
             std::string& physics_table,
             std::string&  physics_dir,
             GlastMS::MultipleScatteringFactory& msfactory,
             GlastMS::EnergyLossFactory& eLossFactory,
			 IG4GeometrySvc*);

    virtual ~RunManager();

  public: 
    /// This is the method to be invoked to start the simulation
    virtual void BeamOn();

    /// Set up the G4 stuff
	virtual void Initialize();

    /// Not used
	virtual void AbortRun(G4bool softAbort = true);

    /// This method handles aborting an event for certain specific circumstances
    virtual void AbortEvent();
  
    /// This method return the number of trajectories stored in the currentEvent
    unsigned int getNumberOfTrajectories();

    /// This method return the vector of points of a trajectory
    int getTrajectoryCharge(unsigned int);

    /// This method return the vector of points of a trajectory
    std::auto_ptr<std::vector<Hep3Vector> > getTrajectoryPoints(unsigned int);

    /// This method return the id (an integer) of the track associated to a
    /// trajectory
    int getTrajectoryTrackId(unsigned int i);

  protected: 
    friend class G4Generator;

    ///  These three protected methods are invoked from Initialize() method
    virtual void InitializeGeometry();
    virtual void InitializePhysics();
    virtual void InitializeCutOff();

    /// Confirm that G4 is in the right state to start an Event
    virtual G4bool ConfirmBeamOnCondition();
    /// Init and termination of a RUN .. to be eliminated?
    virtual void RunInitialization();
    virtual void RunTermination();

	// Method to generate the event
    virtual G4Event* GenerateEvent(G4int i_event);

  public: 

    void DumpRegion(G4String rname) const;
    // Dump information of a region.

    void DumpRegion(G4Region* region=0) const;
    // Dump information of a region.
    // If the pointer is NULL, all regions are shown.

  protected:
    G4RunManagerKernel * kernel;
    G4EventManager * eventManager;

    G4VModularPhysicsList * physicsList;
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

    G4Run*   currentRun;
    G4Event* currentEvent;

    G4int    storeRandomNumberStatus;
    G4String randomNumberStatusDir;

    /// This is a dummy session to be used to silent G4
    UIsession* session;

    /// Range cutoff values for regions
    G4double defaultCut;
    G4double TkrCutValue;
    G4double CalCutValue;

    /// Keep the geometry service for now
    IG4GeometrySvc* m_gsv;

  public:
    virtual void StoreRandomNumberStatus(G4int eventID=-1);
    virtual void RestoreRandomNumberStatus(G4String fileN);

  public: // with description
    //  These methods store respective user initialization and action classes.
    inline const G4VModularPhysicsList* GetUserPhysicsList() const
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
    //  This method must be invoked (or equivalent UI command can be used)
    // in case the user changes his/her detector geometry after
    // Initialize() metho has been invoked. Then, at the begining of the next BeamOn(),
    // all necessary re-initialization will be done.

  public:
    inline void SetVerboseLevel(G4int vl)
    { verboseLevel = vl; 
      kernel->SetVerboseLevel(vl); }
    inline G4int GetVerboseLevel() const
    { return verboseLevel; }

    inline void SetGeometryToBeOptimized(G4bool vl)
    { 
      if(geometryToBeOptimized != vl)
      {
        geometryToBeOptimized = vl;
        geometryNeedsToBeClosed = true;
        kernel->SetGeometryToBeOptimized(vl);
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

