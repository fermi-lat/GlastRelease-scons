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
#include "Event/MonteCarlo/McPositionHit.h"

#include "G4Track.hh"

// This is the singleton static pointer
McParticleManager* McParticleManager::m_pointer = 0;

void McParticleManager::initialize(IDataProviderSvc* esv, IGlastDetSvc* gsvc)
{
	m_esv = esv;
	m_glastDetSvc = gsvc;
}

void McParticleManager::clearRelTables()
{
    // It can happen that a set of relations were stored but never transferred
    // to the TDS because the particle was then rejected. This method will delete
    // any that were not transferred. 

    //**************************************************************************
    //
    // What needs to happen here is that when we delete these relations we need 
    // to also modify the pointers to the parent particle since it is getting 
    // deleted as well... try to set the pointer to the parent of this particle
    //
    //**************************************************************************

    // Start with McParticle to McPositionHit relations
    Event::McPartToPosHitTabList::iterator posHitIter;
    for(posHitIter = m_partToPosHit.begin(); posHitIter != m_partToPosHit.end(); posHitIter++)
    {
        // Clean up the McParticle pointer before deleting relation
        Event::McPartToPosHitRel* rel    = *posHitIter;
        Event::McPositionHit*     posHit = rel->getSecond();
        Event::McParticle*        mcPart = rel->getFirst();

        if (mcPart != mcPart->getMother()) posHit->setMcParticle(mcPart->getMother());
        else                               posHit->setMcParticle(0);

        delete rel;
    }
    m_partToPosHit.clear();

    // same for McParticle to McIntegratingHit
    Event::McPartToIntHitTabList::iterator intHitIter;
    for(intHitIter = m_partToIntHit.begin(); intHitIter != m_partToIntHit.end(); intHitIter++)
    {
        delete *intHitIter;
    }
    m_partToIntHit.clear();

    return;
}


Event::McParticle* McParticleManager::getMcParticle(unsigned int id)
{
  Event::McParticle* particle = 0;

  std::map<unsigned int, Event::McParticle*>::iterator idPartIter = m_particles.find(id);

  if (idPartIter != m_particles.end()) particle = idPartIter->second;

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

    // Should already be clear but just in case
    clearRelTables();
}

void McParticleManager::addMcParticle(unsigned int id, Event::McParticle *particle){
  // Purpose and Method: with that method we add a new McParticle to the map of
  // particles, using as the index an unsigned integer (that is the Geant4 id)

  // We add the particle to the map
  m_particles[id]=particle;
  // We set this particle as the last one added to the map
  m_lastParticle = particle;
  // New particle, start with new relational table list
  clearRelTables();
}

void McParticleManager::addMcRelation(Event::McParticle *particle, Event::McPositionHit* posHit)
{
    Event::McPartToPosHitRel* partToPosHit = new Event::McPartToPosHitRel(particle, posHit);

    m_partToPosHit.push_back(partToPosHit);
}

void McParticleManager::addMcRelation(Event::McParticle *particle, Event::McIntegratingHit* intHit)
{
    Event::McPartToIntHitRel* partToIntHit = new Event::McPartToIntHitRel(particle, intHit);

    m_partToIntHit.push_back(partToIntHit);
}

McParticleManager* McParticleManager::getPointer()
{
  // Purpose and Method: standard singleton method to retrive the unique pointer
  if(m_pointer == 0)
    m_pointer = new McParticleManager;
  return m_pointer;
}

void McParticleManager::clear()
{
    m_particles.clear();
    m_partToPosHit.clear();
    m_partToIntHit.clear();

    // if running FluxAlg, collection will already have parent 
    m_particleColTDS = SmartDataPtr<Event::McParticleCol>(m_esv, EventModel::MC::McParticleCol);

    if( m_particleColTDS == 0) 
    {
        // create the TDS stuff
        m_particleColTDS = new Event::McParticleCol;
        m_esv->registerObject(EventModel::MC::McParticleCol, m_particleColTDS);
    }

    // Try to retrieve the McTrajectoryPoint to McPositionHit relational table
    m_partToPosHitTDS = SmartDataPtr<Event::McPartToPosHitTabList>(m_esv, "/Event/MC/McPartToPosHit");

    if (m_partToPosHitTDS == 0)
    {
        // Create new one and give to TDS
        m_partToPosHitTDS = new Event::McPartToPosHitTabList();
        m_esv->registerObject("/Event/MC/McPartToPosHit", m_partToPosHitTDS);
    }

    // Try to retrieve the McTrajectoryPoint to McIntegratingHit relational table
    m_partToIntHitTDS = SmartDataPtr<Event::McPartToIntHitTabList>(m_esv, "/Event/MC/McPartToIntHit");

    if (m_partToIntHitTDS == 0)
    {
        // Create new one and give to TDS
        m_partToIntHitTDS = new Event::McPartToIntHitTabList();
        m_esv->registerObject("/Event/MC/McPartToIntHit", m_partToIntHitTDS);
    }
}

void McParticleManager::save()
{
    // Purpose and Method: this one save the McParticle hierarchy in the TDS
    // TDS Outputs: the McParticle hierarchy is saved in the
    // /Event/MC/McParticleCol folder of the TDS
  
    // fill the McParticleCol with McParticles
    std::map <unsigned int, Event::McParticle*>::iterator it;
}

bool McParticleManager::makeMcParticle(const G4Track* aTrack)
{
    bool save = true;

    // Search for particle, parentage and process info
    Event::McParticle* parent   = getMcParticle(aTrack->GetParentID()); // Non-zero unless pruning
    Event::McParticle* particle = getMcParticle(aTrack->GetTrackID());  // Is probably zero
  
    // Determine if we should create and save the McParticle info if in pruning mode
    if ( m_mode == PruneMode::MINIMAL_TREE )
    {
        // In full prune mode, we create and save McParticles if
        // 1) This is the primary track
        // 2) This is the direct descendant of the primary track
        if      (particle != 0 && (particle->statusFlags() & Event::McParticle::PRIMARY) != 0) save = true;
        else if (parent   != 0 && (parent->statusFlags()   & Event::McParticle::PRIMARY) != 0) 
        {
            // If a charged track then check to see if it created any hits 
//            if (aTrack->GetDynamicParticle()->GetCharge() != 0.)
//            {
//                if (aTrack->GetKineticEnergy() < m_cutOffE) save = false;
//            }
        }
        else save = false;
    }
    // Save all McParticles by generation (if above cutoff energy)
    else if (m_mode == PruneMode::N_GENERATIONS && aTrack->GetTrackID() > 1)
    {
        save = false;

        if (parent != 0 && aTrack->GetKineticEnergy() > m_cutOffE)
        {
            int genNum = 1;
            SmartRef<Event::McParticle> mcPart = parent->getMother();

            while(mcPart != mcPart->getMother())
            {
                mcPart = mcPart->getMother();
                genNum++;
            }

            if (genNum < m_maxGen) save = true;
        }
    }

    return save;
}

bool McParticleManager::keepMcParticle(const G4Track* aTrack)
{
    // Retrieve pointer to the particle
    Event::McParticle* particle = getMcParticle(aTrack->GetTrackID());  // Is probably zero

    // Ok, if no particle then nothing to do (but shouldn't be called in that case)
    if (!particle) return false;

    // Ok, if a particle then presume we will save it
    bool save = true;

    // If we are in minimal pruning mode, check that this is an "important" McParticle
    if ( m_mode == PruneMode::MINIMAL_TREE )
    {
        // McParticle must have a parent which is the primary particle
        const Event::McParticle* parent = particle->getMother();

        // If the primary we keep it, otherwise check that the secondaries are "important"
        if (parent != 0 && (parent->statusFlags()   & Event::McParticle::PRIMARY) != 0
                        && (particle->statusFlags() & Event::McParticle::PRIMARY) == 0 ) 
        {
            // If a charged track then check to see if it created any hits 
            if (aTrack->GetDynamicParticle()->GetCharge() != 0.)
            {
                if (m_partToPosHit.empty() && m_partToIntHit.empty()) save = false;
            }
        }
    }
/*
    // Save all McParticles by generation (if above cutoff energy)
    else if (m_mode == PruneMode::N_GENERATIONS && aTrack->GetTrackID() > 1)
    {
        save = true;
    }
    // Check for PruneCal mode here 
    else if (m_mode == PruneMode::PRUNE_CAL)
    {
        int towersId;
        int calId;
        m_glastDetSvc->getNumericConstByName("eLATTowers", &towersId);
        m_glastDetSvc->getNumericConstByName("eTowerCAL", &calId);

        if ((particle->getInitialId().size() > 3) && 
	        (particle->getInitialId()[0] == towersId) && 
	        (particle->getInitialId()[3] == calId) &&
	        (particle->getFinalId().size() > 3) &&
	        (particle->getFinalId()[0] == 0) && 
	        (particle->getFinalId()[3] == 0) && 
	       !(particle->statusFlags()&Event::McParticle::POSHIT) && 
	       !(particle->statusFlags()&Event::McParticle::PRIMARY))
        {
            SmartRef<Event::McParticle> mother = &particle->mother();

            mother->removeDaughter(particle);

            save = false;
        }
    }
*/
    return save;
}

void McParticleManager::saveMcParticle()
{
    // Safety just to be sure the particle exists
    if (m_lastParticle)
    {
        // Store McParticle in TDS
        m_particleColTDS->push_back(m_lastParticle);
        m_lastParticle = 0;

        // Move associated relations to the TDS
        Event::McPartToPosHitTabList::iterator posHitIter;
        for(posHitIter = m_partToPosHit.begin(); posHitIter != m_partToPosHit.end(); posHitIter++)
        {
            m_partToPosHitTDS->push_back(*posHitIter);
        }
        m_partToPosHit.clear();

        // Move associated relations to the TDS
        Event::McPartToIntHitTabList::iterator intHitIter;
        for(intHitIter = m_partToIntHit.begin(); intHitIter != m_partToIntHit.end(); intHitIter++)
        {
            m_partToIntHitTDS->push_back(*intHitIter);
        }
        m_partToIntHit.clear();
    }

    return;
}

void McParticleManager::dropMcParticle(const G4Track* aTrack)
{
    // Make sure the particle we want to kill exists
    if (m_lastParticle)
    {
        // Remove this McParticle from map
        std::map<unsigned int, Event::McParticle*>::iterator idPartIter = m_particles.find(aTrack->GetTrackID());
        if (idPartIter != m_particles.end())
        {
            Event::McParticle* mcPart = idPartIter->second;

            // Remove from the daughter list
            const Event::McParticle* mother = mcPart->getMother();
            const_cast<Event::McParticle*>(mother)->removeDaughter(mcPart);

            if (m_lastParticle != mcPart) 
            {
                int j = 0;
            }
            delete mcPart;

            m_particles.erase(idPartIter);
        }

        m_lastParticle = 0;
        clearRelTables();
    }

    return;
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






