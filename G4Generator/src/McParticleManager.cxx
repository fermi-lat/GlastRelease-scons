#include "McParticleManager.h"

// Gaudi
#include "GaudiKernel/SmartDataPtr.h"


McParticleManager* McParticleManager::m_pointer = 0;

McParticleManager* McParticleManager::getPointer()
{
  if(m_pointer == 0)
    m_pointer = new McParticleManager;
  return m_pointer;
}

void McParticleManager::save()
{
  // create the TDS stuff
  mc::McParticleCol* pcol = new mc::McParticleCol;
  m_esv->registerObject("/Event/MC/McParticleCol", pcol);

  // fill the McParticleCol with McParticles
  std::map <unsigned int, mc::McParticle*>::iterator it;

  for(it=m_particles.begin();it != m_particles.end(); it++)
    pcol->push_back(it->second);
}








