#include "McParticleManager.h"

// Gaudi
#include "GaudiKernel/SmartDataPtr.h"


McParticleManager* McParticleManager::m_pointer = 0;

void McParticleManager::addMcParticle(unsigned int id, 
                                      mc::McParticle *particle){
  m_particles[id]=particle; 
  m_lastParticle = particle;
}

McParticleManager* McParticleManager::getPointer()
{
  if(m_pointer == 0)
    m_pointer = new McParticleManager;
  return m_pointer;
}

void McParticleManager::save()
{
    // if running FluxAlg, collection will already have parent 
    mc::McParticleCol*  pcol=  SmartDataPtr<mc::McParticleCol>(m_esv, "/Event/MC/McParticleCol");

    if( pcol==0) {
        // create the TDS stuff
        pcol = new mc::McParticleCol;
        m_esv->registerObject("/Event/MC/McParticleCol", pcol);
    }

    // fill the McParticleCol with McParticles
  std::map <unsigned int, mc::McParticle*>::iterator it;

  for(it=m_particles.begin();it != m_particles.end(); it++)
    pcol->push_back(it->second);
}








