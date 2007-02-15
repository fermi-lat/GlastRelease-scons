#ifndef McTrajectoryManager_H
#define McTrajectoryManager_H

#include <algorithm>

//#include "GaudiKernel/Algorithm.h"
//#include "GaudiKernel/MsgStream.h"

#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McTrajectory.h"
#include "Event/MonteCarlo/McRelTableDefs.h"
#include "Event/RelTable/RelTable.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include <map>

class IGlastDetSvc;
class G4Track;

/** 
 *  @class McTrajectoryManager
 *
 *  @brief A singleton manager for McParticle hierarchy
 *
 *  The McParticle tree is built in two different places: in the G4Generator,
 *  where it is initialized with the primary particle, and in the TrackingAction
 *  where new particles are added to the tree with relevant parentage
 *  information. Moreover the McParticle hierarchy must be accessed during hit
 *  filling of the TDS by the detector managers. For these reasons this
 *  McTrajectoryManager class has been implemented: being a singleton only one
 *  instance is created and can be accessed everywhere it is needed.
 *
 *  @author R.Giannitrapani
 *
 * $Header$
 */
class McTrajectoryManager {
 public:
  /// The static pointer retrival method of the singleton
  static McTrajectoryManager* getPointer();

  /// This method adds an McTrajectory collection to the map with
  /// @param id The Geant4 id (a progressive integer)
  /// @param particle The pointer to the McParticle to be added
  void addMcTrajectory(Event::McTrajectoryCol *pcol);

  /// This method add an McParticle to the map with id as an index
  /// @param id The Geant4 id (a progressive integer)
  /// @param particle The pointer to the McParticle to be added
  /// @param optional pointer to McParticle for relational table
  void addMcTrajectory(unsigned int id, Event::McTrajectory *trajectory, Event::McParticle* mcPart=0);

  /// This method sets the McPositionHit This method is to be called 
  /// by the PosDetectorManager during stepping
  /// @param mcPorIHit pointer to the hit in question
  void setMcPosHit(const Event::McPositionHit* mcPosHit) {m_mcPosHit = mcPosHit;}

  /// This method returns the ContainedObject pointer set above
  const Event::McPositionHit* getMcPosHit() const {return m_mcPosHit;}

  /// This method sets the McTrajectoryPoint to McPositionHit relation. Is called
  /// by the SteppingAction called at the end of a step. 
  /// @param mcPorIHit pointer to the hit in question
  void addMcPosRel(Event::McPointToPosHitRel* mcPosHitRel) {m_pointToPosHit.push_back(mcPosHitRel);}

  /// This method sets the McIntegratingHit. This method is to be called 
  /// by the IntDetectorManager during stepping
  /// @param mcPorIHit pointer to the hit in question
  void setMcIntHit(const Event::McIntegratingHit* mcIntHit) {m_mcIntHit = mcIntHit;}

  /// This method returns the McIntegratingHit pointer set above
  const Event::McIntegratingHit* getMcIntHit() const {return m_mcIntHit;}

  /// This method sets the McTrajectoryPoint to McIntegratingHit relation. Is called
  /// by the SteppingAction called at the end of a step. 
  /// @param mcPorIHit pointer to the hit in question
  void addMcIntRel(Event::McPointToIntHitRel* mcIntHitRel) {m_pointToIntHit.push_back(mcIntHitRel);}

  /// Retrive an McParticle giving an id
  /// @param id the Geant4 id (a progressive integer)
  /// @return The McParticle corresponding to the id
  Event::McTrajectory* getMcTrajectory(unsigned int id);

  /// This method is called at the end of tracking a single particle and
  /// will move the McTrajectory and relational table information
  /// cached for this track into the TDS 
  void saveMcTrajectory();

  /// This method will drop the McTrajectory and relational table information
  /// cached for this track into the TDS 
  void dropMcTrajectory(const G4Track* track);

  /// initialize the class, passing the event selector pointer for TDS
  /// operations
  /// @param esv The pointer to the DataProviderSvc service 
  /// @param gsvc The pointer to the GlastDetSvc service
  void initialize(IDataProviderSvc* esv, IGlastDetSvc* gsvc);
  
  /// clear the hierarchy of McTrajectory
  void clear();

  /// Save the McTrajectory hierarchy in the TDS
  void save();

  /// Return the number of McTrajectory's saved
  unsigned int size(){return m_trajectories.size();};

  /// For pruning, remove a given trajectory
  void removeTrajectory(unsigned int id);

 private:
  /// The constructor is private since this is a singleton
  McTrajectoryManager() : 
      m_mcTrajectory(0), m_trajectoryColTDS(0), m_pointToPosHitTDS(0), m_pointToIntHitTDS(0), m_mcPosHit(0), m_mcIntHit(0)
      {}   

  /// private method to keep the local relational table lists clear
  void clearRelTables();

  /// The static pointer of the singleton
  static McTrajectoryManager* m_pointer;
  
  /// The map of McTrajectory pointers indicized by g4 ids 
  std::map <unsigned int, Event::McTrajectory*> m_trajectories;

  /// Local copy of the actual McTrajectory
  Event::McTrajectory* m_mcTrajectory;

  /// Local copy of McParticle to McTrajectory relation (there will be one per track)
  Event::McPartToTrajectoryRel* m_partToTrajectory;

  /// Local list of McPoint to McPositionHits
  std::list<Event::McPointToPosHitRel*> m_pointToPosHit;

  /// Local list of McPoint to McPositionHits
  std::list<Event::McPointToIntHitRel*> m_pointToIntHit;

  /// TDS Container for the McTrajectory objects
  Event::McTrajectoryCol* m_trajectoryColTDS;

  /// TDS Relational table between McParticles and McTrajectorys
  Event::McPartToTrajectoryTabList* m_partToTrajectoryTDS;

  /// TDS Relational table betweeen trajectory points and McPositionHits
  Event::McPointToPosHitTabList* m_pointToPosHitTDS;

  /// TDS Relational table between trajectory points and McIntegratingHits
  Event::McPointToIntHitTabList* m_pointToIntHitTDS;

  /// Pointer to either an McIntegratingHit or McPositionHit
  const Event::McPositionHit*    m_mcPosHit;
  const Event::McIntegratingHit* m_mcIntHit;

  /// The pointer to the IdataProviderSvc
  IDataProviderSvc* m_esv;

  /// The GlastDetSvc is used to retrive some constants in the pruneCal method 
  IGlastDetSvc* m_glastDetSvc;
};
#endif //McTrajectoryManager_H
