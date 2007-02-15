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

#include "McTrajectoryManager.h"

// Gaudi
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRef.h"
#include "GaudiKernel/SmartRefVector.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "Event/TopLevel/EventModel.h"

// This is the singleton static pointer
McTrajectoryManager* McTrajectoryManager::m_pointer = 0;

void McTrajectoryManager::initialize(IDataProviderSvc* esv, IGlastDetSvc* gsvc)
{
    m_mcPosHit            = 0;
    m_mcIntHit            = 0;
	m_esv                 = esv;
	m_glastDetSvc         = gsvc;
    m_pointToPosHitTDS    = 0;
    m_pointToIntHitTDS    = 0;
    m_partToTrajectoryTDS = 0;
    m_mcTrajectory        = 0;
    m_partToTrajectory    = 0;

    m_trajectories.clear();
    m_pointToPosHit.clear();
    m_pointToIntHit.clear();
}

void McTrajectoryManager::clearRelTables()
{
    // It can happen that a set of relations were stored but never transferred
    // to the TDS because the particle was then rejected. This method will delete
    // any that were not transferred. 

    // Not a table but a releation
    if (m_partToTrajectory) delete m_partToTrajectory;
    m_partToTrajectory = 0;

    // Start with McTrajectoryPoint to McPositionHit relations
    std::list<Event::McPointToPosHitRel*>::iterator posHitIter;
    for(posHitIter = m_pointToPosHit.begin(); posHitIter != m_pointToPosHit.end(); posHitIter++)
    {
        delete *posHitIter;
    }
    m_pointToPosHit.clear();

    // Start with McTrajectoryPoint to McIntegratingHit relations
    std::list<Event::McPointToIntHitRel*>::iterator intHitIter;
    for(intHitIter = m_pointToIntHit.begin(); intHitIter != m_pointToIntHit.end(); intHitIter++)
    {
        delete *intHitIter;
    }
    m_pointToIntHit.clear();

    return;
}

void McTrajectoryManager::clear()
{
    m_mcPosHit     = 0;
    m_mcIntHit     = 0;
    m_mcTrajectory = 0;
    m_trajectories.clear(); 

    clearRelTables();

    // Try to get the McTrajectory collection from the TDS 
    m_trajectoryColTDS = SmartDataPtr<Event::McTrajectoryCol>(m_esv, EventModel::MC::McTrajectoryCol);

    if( m_trajectoryColTDS == 0) 
    {
        // create the TDS stuff
        m_trajectoryColTDS = new Event::McTrajectoryCol;
        m_esv->registerObject(EventModel::MC::McTrajectoryCol, m_trajectoryColTDS);
    }

    // Try to retrieve the McParticle to McTrajectory relational table
    m_partToTrajectoryTDS = SmartDataPtr<Event::McPartToTrajectoryTabList>(m_esv, "/Event/MC/McPartToTrajectory");

    if (m_partToTrajectoryTDS == 0)
    {
        // Create new one and give to TDS
        m_partToTrajectoryTDS = new Event::McPartToTrajectoryTabList();
        m_esv->registerObject("/Event/MC/McPartToTrajectory", m_partToTrajectoryTDS);
    }

    // Try to retrieve the McTrajectoryPoint to McPositionHit relational table
    m_pointToPosHitTDS = SmartDataPtr<Event::McPointToPosHitTabList>(m_esv, "/Event/MC/McPointToPosHit");

    if (m_pointToPosHitTDS == 0)
    {
        // Create new one and give to TDS
        m_pointToPosHitTDS = new Event::McPointToPosHitTabList();
        m_esv->registerObject("/Event/MC/McPointToPosHit", m_pointToPosHitTDS);
    }

    // Try to retrieve the McTrajectoryPoint to McIntegratingHit relational table
    m_pointToIntHitTDS = SmartDataPtr<Event::McPointToIntHitTabList>(m_esv, "/Event/MC/McPointToIntHit");

    if (m_pointToIntHitTDS == 0)
    {
        // Create new one and give to TDS
        m_pointToIntHitTDS = new Event::McPointToIntHitTabList();
        m_esv->registerObject("/Event/MC/McPointToIntHit", m_pointToIntHitTDS);
    }
}


Event::McTrajectory* McTrajectoryManager::getMcTrajectory(unsigned int id)
{
  Event::McTrajectory* trajectory = 0;

  std::map <unsigned int, Event::McTrajectory*>::iterator trajIter = m_trajectories.find(id);

  if (trajIter != m_trajectories.end()) trajectory = trajIter->second;

  return trajectory;
}

void McTrajectoryManager::addMcTrajectory(Event::McTrajectoryCol *pcol)
{
    // Purpose and Method: with that method we add a new McTrajectory to the map of
    // particles, using as the index an unsigned integer (that is the Geant4 id)
    // McParticle counter...
    int idx = 0;

    Event::McTrajectoryCol::iterator pColIter;

    for(pColIter = pcol->begin(); pColIter != pcol->end(); pColIter++)
    {
        Event::McTrajectory* mcTrajectory = *pColIter;
        
        addMcTrajectory(idx++, mcTrajectory);
    }
}

void McTrajectoryManager::addMcTrajectory(unsigned int id, Event::McTrajectory *trajectory, Event::McParticle* mcPart){
  // Purpose and Method: with that method we add a new McParticle to the map of
  // particles, using as the index an unsigned integer (that is the Geant4 id)

  // Add the trajectory to the collection (which will be THE owner)
  //m_trajectoryCol->push_back(trajectory);
  m_mcTrajectory = trajectory;

  // Add to relational table if McParticle pointer supplied
  if (mcPart)
  {
      m_partToTrajectory = new Event::McPartToTrajectoryRel(mcPart, trajectory);
  }
  else m_partToTrajectory = 0;

  // We add the trajectory to the map
  m_trajectories[id]=trajectory; 

  // Clear the pointers caching hit pointers (since this is new particle)
  m_mcPosHit = 0;
  m_mcIntHit = 0;

}
void McTrajectoryManager::saveMcTrajectory()
{
    // Make sure we have something to store
    if (m_mcTrajectory)
    {
        // Move the trajectory itself into the TDS
        m_trajectoryColTDS->push_back(m_mcTrajectory);
        m_mcTrajectory = 0;

        // The McParticle to McTrajectory relation
        m_partToTrajectoryTDS->push_back(m_partToTrajectory);
        m_partToTrajectory = 0;

        // Move associated relations to the TDS
        std::list<Event::McPointToPosHitRel*>::iterator posHitIter;
        for(posHitIter = m_pointToPosHit.begin(); posHitIter != m_pointToPosHit.end(); posHitIter++)
        {
            m_pointToPosHitTDS->push_back(*posHitIter);
        }
        m_pointToPosHit.clear();

        // Move associated relations to the TDS
        std::list<Event::McPointToIntHitRel*>::iterator intHitIter;
        for(intHitIter = m_pointToIntHit.begin(); intHitIter != m_pointToIntHit.end(); intHitIter++)
        {
            m_pointToIntHitTDS->push_back(*intHitIter);
        }
        m_pointToIntHit.clear();
    }
}

void McTrajectoryManager::dropMcTrajectory(const G4Track* track)
{
    // Make sure we have something to get rid of
    if (m_mcTrajectory)
    {
        clearRelTables();

        delete m_mcTrajectory;

        m_mcTrajectory = 0;
    }
}

McTrajectoryManager* McTrajectoryManager::getPointer()
{
  // Purpose and Method: standard singleton method to retrive the unique pointer
  if(m_pointer == 0)
    m_pointer = new McTrajectoryManager;
  return m_pointer;
}

void McTrajectoryManager::save()
{
    // Purpose and Method: this one save the McParticle hierarchy in the TDS
    // TDS Outputs: the McParticle hierarchy is saved in the
    // /Event/MC/McParticleCol folder of the TDS

    // Move along, nothing to see here!! 
    return;
}

void McTrajectoryManager::removeTrajectory(unsigned int id)
{
    m_trajectories.erase(id);
}




