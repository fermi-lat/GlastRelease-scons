#ifndef MCPARTICLEMANAGER_H
#define MCPARTICLEMANAGER_H

#include <algorithm>

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McRelTableDefs.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include <map>

class IGlastDetSvc;
class G4Track;

/** 
 *  @class McParticleManager
 *
 *  @brief A singleton manager for McParticle hierarchy
 *
 *  The McParticle tree is built in two different places: in the G4Generator,
 *  where it is initialized with the primary particle, and in the TrackingAction
 *  where new particles are added to the tree with relevant parentage
 *  information. Moreover the McParticle hierarchy must be accessed during hit
 *  filling of the TDS by the detector managers. For these reasons this
 *  McParticleManager class has been implemented: being a singleton only one
 *  instance is created and can be accessed everywhere it is needed.
 *
 *  @author R.Giannitrapani
 *
 * $Header$
 */
class McParticleManager {
 public:
  /// Enumerate the available McParticle pruning modes
  enum PruneMode {
                   FULL_MC_TREE,    // No pruning
                   MINIMAL_TREE,    // Keep only primary and immediate daughters
                   N_GENERATIONS,   // Keep up to N generations
                   PRUNE_CAL        // Prune calorimeter
  };
  /// The static pointer retrival method of the singleton
  static McParticleManager* getPointer();

  /// This method add an McParticle to the map with id as an index
  /// @param id The Geant4 id (a progressive integer)
  /// @param particle The pointer to the McParticle to be added
  void addMcParticle(Event::McParticleCol *pcol);

  /// This method add an McParticle to the map with id as an index
  /// @param id The Geant4 id (a progressive integer)
  /// @param particle The pointer to the McParticle to be added
  void addMcParticle(unsigned int id, Event::McParticle *particle);

  /// This method adds an McParticle <--> McPositionHit relation 
  /// to the relations list
  /// @param particle - pointer to the McParticle
  /// @param posHit - pointer to the McPositionHit
  void addMcRelation(Event::McParticle *particle, Event::McPositionHit* posHit);

  /// This method adds an McParticle <--> McIntegratingHit relation 
  /// to the relations list
  /// @param particle - pointer to the McParticle
  /// @param intHit - pointer to the McIntegratingHit
  void addMcRelation(Event::McParticle *particle, Event::McIntegratingHit* intHit);

  /// Retrive an McParticle giving an id
  /// @param id the Geant4 id (a progressive integer)
  /// @return The McParticle corresponding to the id
  Event::McParticle* getMcParticle(unsigned int id);

  /// Retrive the last particle added to the map
  /// @return The last McParticle added to the map
  Event::McParticle* getLastParticle()const{return m_lastParticle;};

  /// Retrive the last particle added to the map
  /// @return The last McParticle added to the map
  //Event::McParticle* setLastParticle(Event::McParticle* particle) {m_lastParticle = particle;};

  /// initialize the class, passing the event selector pointer for TDS
  /// operations
  /// @param esv The pointer to the DataProviderSvc service 
  /// @param gsvc The pointer to the GlastDetSvc service
  void initialize(IDataProviderSvc* esv, IGlastDetSvc* gsvc);

  /// Retrive the actual origin particle 
  /// @return The actual origin McParticle 
  Event::McParticle* getOriginParticle()const{return m_currentOrigin;};

  /// Retrive the actual origin particle id 
  /// @return The actual origin McParticle id
  Event::McParticle::StdHepId getOriginParticleId(){return m_currentOriginId;};

  /// Set the actual origin particle 
  void setOriginParticle(Event::McParticle* origin){m_currentOrigin = origin;};

  /// Set the actual origin particle id 
  void setOriginParticleId(Event::McParticle::StdHepId id){m_currentOriginId = id;};

  /// Set and get methods for the mode of the McParticle tree saving
  void setMode(PruneMode m) {m_mode = m;};
  PruneMode getMode() {return m_mode;};

  /// Set and get methods for the cutoff energy (various pruning modes)
  void setCutOffEnergy(float cutOffE) {m_cutOffE = cutOffE;}
  float getCutOffEnergy() {return m_cutOffE;}

  /// Set and get methods for the number of generations (N_GENERATIONS mode)
  void setNumGenerations(int nGen) {m_maxGen = nGen;}
  float getNumGenerations() {return m_maxGen;}
  
  /// clear the hierarchy of McParticle
  void clear();

  /// Save the McParticle hierarchy in the TDS
  void save();

  /// Return the number of McParticle saved
  unsigned int size(){return m_particles.size();};

  /// At beginning of a track, decide whether to make an McParticle
  /// This depends on pruning mode we are running
  bool makeMcParticle(const G4Track* track);

  /// At end of a track, decide whether to keep the current McParticle
  /// This depends on pruning mode we are running
  bool keepMcParticle(const G4Track* track);

  /// This method is called at the end of tracking a single particle and
  /// will move all McParticle and relational table information
  /// cached for this track into the TDS 
  void saveMcParticle();

  /// This method will drop all McParticle and relational table information
  /// cached for this track into the TDS 
  void dropMcParticle(const G4Track* track);

  /// Prune the tree in the CAL part
  void pruneCal();

 private:
  /// The constructor is private since this is a singleton
     McParticleManager() : 
        m_lastParticle(0),m_currentOrigin(0),m_mode(FULL_MC_TREE),m_cutOffE(0.),m_maxGen(3) 
     {m_particles.clear(); m_partToPosHit.clear(); m_partToIntHit.clear();}; 

  /// private method to keep the local relational table lists clear
  void clearRelTables();

  /// The static pointer of the singleton
  static McParticleManager* m_pointer;
  
  /// The map of McParticle pointers indicized by g4 ids 
  std::map <unsigned int, Event::McParticle*> m_particles;

  /// Local list of McParticle to McPositionHit relations
  std::list<Event::McPartToPosHitRel*> m_partToPosHit;

  /// TDS List of relations between McParticles and McIntegratingHits
  std::list<Event::McPartToIntHitRel*> m_partToIntHit;

  /// Pointer to the TDS collection of McParticles 
  /// These are only the ones we keep
  Event::McParticleCol* m_particleColTDS;

  /// TDS List of relations between McParticles and McPositionHits
  Event::McPartToPosHitTabList* m_partToPosHitTDS;

  /// TDS List of relations between McParticles and McIntegratingHits
  Event::McPartToIntHitTabList* m_partToIntHitTDS;

  /// The pointer to the IdataProviderSvc
  IDataProviderSvc* m_esv;

  /// A pointer to the last particle added to the map
  Event::McParticle* m_lastParticle;

  /// A pointer to the last ancestor (e+ or e- in the case of a primary gamma,
  /// the primary itself otherwise)
  Event::McParticle* m_currentOrigin;

  /// The id of the last ancestor (e+ or e- in the case of a primary gamma,
  /// the primary itself otherwise)
  Event::McParticle::StdHepId m_currentOriginId;

  /// A flag for the mode of generation of the McParticle tree; true means the
  /// full modes false means the minimal one
  PruneMode m_mode;

  /// Lower cut off energy for pruning modes
  float m_cutOffE;

  /// Maximum generations to keep
  int m_maxGen;

  /// The GlastDetSvc is used to retrive some constants in the pruneCal method 
  IGlastDetSvc* m_glastDetSvc;
};
#endif //McParticleManager_H
