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
#include "GaudiKernel/SmartRef.h"
#include "GaudiKernel/SmartRefVector.h"

#include "Event/TopLevel/EventModel.h"

// This is the singleton static pointer
McParticleManager* McParticleManager::m_pointer = 0;

Event::McParticle* McParticleManager::getMcParticle(unsigned int id)
{
  Event::McParticle* particle = 0;
  std::map<unsigned int, Event::McParticle*>::iterator it;
  
  for(it=m_particles.begin();it != m_particles.end() ; it++){
    if (it->first == id) 
      {
        particle = it->second;
        break;
      }
  }

  return particle;
}
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
    if(it->second)
      pcol->push_back(it->second);
  }

}

void McParticleManager::pruneCal()
{
  // Purpose and Method: with that method we prune all the McParticle that does
  // not belong to the TKR region (very LAT speicific, we check z<0 for both
  // start and end position of an McParticle) and that does not interact with the
  // TKR itself (producing an McPositionHit).

  SmartRefVector<Event::McParticle>::const_iterator d; 
  std::map <unsigned int, Event::McParticle*>::iterator it;

  for(it=m_particles.begin();it != m_particles.end() ; it++)
    {
      if (((it->second)->initialPosition().z() < 0) && ((it->second)->finalPosition().z()<0)
          && !((it->second)->statusFlags()&Event::McParticle::POSHIT)
          && !((it->second)->statusFlags()&Event::McParticle::PRIMARY))
        {
          for(unsigned int i=0;i<(it->second)->daughterList().size();i++)
            {
              SmartRef<Event::McParticle> part = (it->second)->daughterList()[i];
              SmartRef<Event::McParticle> mother = &(it->second)->mother();
              part->setMother(mother);
              part->addStatusFlag(Event::McParticle::BCKSPL);
              mother->addDaughter(part);
            }

          delete it->second;
          it->second = 0;
        }
      
    }
}






