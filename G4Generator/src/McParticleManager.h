#ifndef MCPARTICLEMANAGER_H
#define MCPARTICLEMANAGER_H

#include <algorithm>

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"

#include "GlastEvent/MonteCarlo/McParticle.h"
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
  void addMcParticle(unsigned int id, mc::McParticle *particle);

  /// Retrive an McParticle giving an id
  /// @param id the Geant4 id (a progressive integer)
  /// @return The McParticle corresponding to the id
  mc::McParticle* getMcParticle(unsigned int id){return m_particles[id];}

  /// Retrive the last particle added to the map
  /// @return The last McParticle added to the map
  mc::McParticle* getLastParticle(){return m_lastParticle;};

  /// initialize the class, passing the event selector pointer for TDS
  /// operations
  /// @param esv The pointer to the DataProviderSvc service 
  void initialize(IDataProviderSvc* esv){m_esv = esv;}
  
  /// clear the hierarchy of McParticle
  void clear(){m_particles.clear();}

  /// Save the McParticle hierarchy in the TDS
  void save();

 private:
  /// The constructor is private since this is a singleton
  McParticleManager():m_lastParticle(0){}; 

  /// The static pointer of the singleton
  static McParticleManager* m_pointer;
  
  /// The map of McParticle pointers indicized by g4 ids 
  std::map <unsigned int, mc::McParticle*> m_particles;

  /// The pointer to the IdataProviderSvc
  IDataProviderSvc* m_esv;

  /// A pointer to the last particle added to the map
  mc::McParticle* m_lastParticle;
};
#endif //McParticleManager_H
