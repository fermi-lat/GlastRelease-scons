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
 *  @author R.Giannitrapani
 */

class McParticleManager {
 public:

  /// The static pointer retrival method of the singleton
  static McParticleManager* getPointer();

  /// This method add an McParticle to the map with id as an index
  void addMcParticle(unsigned int id, mc::McParticle *particle)
    {m_particles[id]=particle; 
     m_lastParticle = particle;};

  /// Retrive an McParticle giving an id
  mc::McParticle* getMcParticle(unsigned int id){return m_particles[id];}

  /// Retrive the last particle added to the map
  mc::McParticle* getLastParticle(){return m_lastParticle;};

  /// initialize the class, passing the event selector pointer for TDS
  /// operations
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

  mc::McParticle* m_lastParticle;
  };

#endif //McParticleManager_H
