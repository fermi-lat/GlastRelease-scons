// File and Version Information:
// $Header$
//
// Description: this utility singleton is used in various classes of G4Generator
// to register new McParticle objects, retrive the actual McParticle (i.e. the
// one that is actually creating hits in the sensitive detectors) and to finally
// save the McParticle hierarchy in the TDS
//      
// Author(s):
//      R.Giannitrapani

#include "McParticleManager.h"

// Gaudi
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/TopLevel/EventModel.h"

// This is the singleton static pointer
McParticleManager* McParticleManager::m_pointer = 0;

void McParticleManager::addMcParticle(unsigned int id, 
                                      Event::McParticle *particle){
  // Purpose and Method: with that method we add a new McParticle to the map of
  // particles, using as the index an unsigned integer (that is the Geant4 id)

  // We add the particle to the map
  m_particles[id]=particle; 
  // We set this particle as the last one added to the map
  m_lastParticle = particle;
}

McParticleManager* McParticleManager::getPointer()
{
  // Purpose and Method: standard singleton method to retrive the unique pointer
  if(m_pointer == 0)
    m_pointer = new McParticleManager;
  return m_pointer;
}

void McParticleManager::save()
{
  // Purpose and Method: this one save the McParticle hierarchy in the TDS
  // TDS Outputs: the McParticle hierarchy is saved in the
  // /Event/MC/McParticleCol folder of the TDS

  // if running FluxAlg, collection will already have parent 
  Event::McParticleCol*  pcol=  
    SmartDataPtr<Event::McParticleCol>(m_esv, EventModel::MC::McParticleCol);

  if( pcol==0) {
    // create the TDS stuff
    pcol = new Event::McParticleCol;
    m_esv->registerObject(EventModel::MC::McParticleCol, pcol);
  }
  
  // fill the McParticleCol with McParticles
  std::map <unsigned int, Event::McParticle*>::iterator it;
  
  for(it=m_particles.begin();it != m_particles.end() ; it++){
    pcol->push_back(it->second);
  }
}








