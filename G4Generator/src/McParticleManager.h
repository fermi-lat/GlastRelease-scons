#ifndef MCPARTICLEMANAGER_H
#define MCPARTICLEMANAGER_H

#include <algorithm>

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/MonteCarlo/McParticle.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include <map>

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
  /// The static pointer retrival method of the singleton
  static McParticleManager* getPointer();

  /// This method add an McParticle to the map with id as an index
  /// @param id The Geant4 id (a progressive integer)
  /// @param particle The pointer to the McParticle to be added
  void addMcParticle(unsigned int id, Event::McParticle *particle);

  /// Retrive an McParticle giving an id
  /// @param id the Geant4 id (a progressive integer)
  /// @return The McParticle corresponding to the id
  Event::McParticle* getMcParticle(unsigned int id);

  /// Retrive the last particle added to the map
  /// @return The last McParticle added to the map
  Event::McParticle* getLastParticle()const{return m_lastParticle;};

  /// initialize the class, passing the event selector pointer for TDS
  /// operations
  /// @param esv The pointer to the DataProviderSvc service 
  void initialize(IDataProviderSvc* esv){m_esv = esv;}

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
  void setMode(bool m){m_mode = m;};
  bool getMode(){return m_mode;};
  
  /// clear the hierarchy of McParticle
  void clear(){m_particles.clear();}

  /// Save the McParticle hierarchy in the TDS
  void save();

  /// Return the number of McParticle saved
  unsigned int size(){return m_particles.size();};

  /// Prune the tree in the CAL part
  void pruneCal();

 private:
  /// The constructor is private since this is a singleton
  McParticleManager():m_lastParticle(0),m_currentOrigin(0),m_mode(1){clear();}; 

  /// The static pointer of the singleton
  static McParticleManager* m_pointer;
  
  /// The map of McParticle pointers indicized by g4 ids 
  std::map <unsigned int, Event::McParticle*> m_particles;

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
  bool m_mode;
};
#endif //McParticleManager_H
