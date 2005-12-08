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

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "Event/TopLevel/EventModel.h"

// This is the singleton static pointer
McParticleManager* McParticleManager::m_pointer = 0;

void McParticleManager::initialize(IDataProviderSvc* esv, IGlastDetSvc* gsvc)
{
	m_esv = esv;
	m_glastDetSvc = gsvc;
}


Event::McParticle* McParticleManager::getMcParticle(unsigned int id)
{
  Event::McParticle* particle = 0;
  
  particle = m_particles[id];

  //std::map<unsigned int, Event::McParticle*>::iterator it;
  //for(it=m_particles.begin();it != m_particles.end() ; it++){
  //  if (it->first == id) 
  //    {
  //      particle = it->second;
  //      break;
  //    }
  //}

  return particle;
}
void McParticleManager::addMcParticle(Event::McParticleCol *pcol)
{
    // Purpose and Method: with that method we add a new McParticle to the map of
    // particles, using as the index an unsigned integer (that is the Geant4 id)
    // McParticle counter...
    int idx = 0;

    // The primary particle
    Event::McParticle* primary = pcol->front();

    // Add the primary to our McParticle <--> id map
    addMcParticle(idx++, primary);

    // If the primary has daughters then we need to add em too!
    // Retrieve the vector of daughter particles
    const SmartRefVector<Event::McParticle>& daughterVec = primary->daughterList();
  
    SmartRefVector<Event::McParticle>::const_iterator dIter = daughterVec.begin();

    int numDaughters = daughterVec.size();

    //for(Event::McParticleCol::iterator colIter = pcol->begin(); colIter != pcol->end(); colIter++)
    for(dIter = daughterVec.begin(); dIter != daughterVec.end(); dIter++)
    {
        const Event::McParticle* mcPart  = *dIter;

        addMcParticle(idx++, const_cast<Event::McParticle*>(mcPart));
    }
}

void McParticleManager::addMcParticle(unsigned int id, Event::McParticle *particle){
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

    if( pcol==0) 
    {
        // create the TDS stuff
        pcol = new Event::McParticleCol;
        m_esv->registerObject(EventModel::MC::McParticleCol, pcol);
    }
  
    // fill the McParticleCol with McParticles
    std::map <unsigned int, Event::McParticle*>::iterator it;
  
    for(it=m_particles.begin();it != m_particles.end() ; it++)
    {
        if(it->second) pcol->push_back(it->second);
    }
}

void McParticleManager::pruneCal()
{
  // Purpose and Method: with that method we prune all the McParticle that 
  // start in the Cal and die in the Cal (to find that we are using the id
  // and some GlastDetSvc constant to check them.

  std::map <unsigned int, Event::McParticle*>::iterator it;

  int towersId;
  int calId;
  m_glastDetSvc->getNumericConstByName("eLATTowers", &towersId);
  m_glastDetSvc->getNumericConstByName("eTowerCAL", &calId);

  for(it=m_particles.begin();it != m_particles.end() ; it++)
    {
      SmartRef<Event::McParticle> mcPart = it->second;

      if ((mcPart->getInitialId().size() > 3) && 
	      (mcPart->getInitialId()[0] == towersId) && 
	      (mcPart->getInitialId()[3] == calId) &&
	      (mcPart->getFinalId().size() > 3) &&
	      (mcPart->getFinalId()[0] == 0) && 
	      (mcPart->getFinalId()[3] == 0) && 
	     !(mcPart->statusFlags()&Event::McParticle::POSHIT) && 
	     !(mcPart->statusFlags()&Event::McParticle::PRIMARY))
        {
          SmartRef<Event::McParticle> mother = &mcPart->mother();

          mother->removeDaughter(mcPart);

          for(unsigned int i=0;i<mcPart->daughterList().size();i++)
            {
              SmartRef<Event::McParticle> part = mcPart->daughterList()[i];
              part->setMother(mother);
              part->addStatusFlag(Event::McParticle::BCKSPL);
              mother->addDaughter(part);
            }

          delete &(*mcPart); // unix expects a real pointer
          it->second = 0;
        }
      
    }
}






