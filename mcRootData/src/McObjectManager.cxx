// File and Version Information:
// $Header$
//
// Description: this utility singleton is used in various classes of G4Generator
// to register new McParticle objects, retrive the actual McParticle (i.e. the
// one that is actually creating hits in the sensitive detectors) and to finally
// save the McParticle hierarchy in the TDS
//      
// Author(s):
//      Am. Crazy

#include "McObjectManager.h"

// This is the singleton static pointer
McObjectManager* McObjectManager::m_pointer = 0;

// Define the pool for McParticles
#define MCPOOLSIZE         1000
#define MCMAXPOOLGROW      50000
#define MCHITPOOLSIZE      1000
#define MCHITMAXPOOLGROW   1000
#define MCPOINTPOOLSIZE    5000
#define MCPOINTMAXPOOLGROW 50000

McObjectManager::McObjectManager() : m_mcPartPool(MCPOOLSIZE), 
                                     m_mcPosHitPool(MCHITPOOLSIZE), 
                                     m_mcIntHitPool(MCHITPOOLSIZE),
                                     m_mcTrajectoryPool(MCPOOLSIZE),
                                     m_mcTrajPointPool(MCPOINTPOOLSIZE)
{
    m_mcPartPoolIdx       = m_mcPartPool.begin();
    m_mcPosHitPoolIdx     = m_mcPosHitPool.begin();
    m_mcIntHitPoolIdx     = m_mcIntHitPool.begin();
    m_mcTrajectoryPoolIdx = m_mcTrajectoryPool.begin();
    m_mcTrajPointPoolIdx  = m_mcTrajPointPool.begin();
}

McObjectManager* McObjectManager::getPointer()
{
  // Purpose and Method: standard singleton method to retrive the unique pointer
  if(m_pointer == 0) m_pointer = new McObjectManager();
  return m_pointer;
}

McParticle* McObjectManager::getNewMcParticle()
{
    McParticle* mcPart = 0;

    // If we have exceeded our pre-allocated list of McParticles then expand
    if (m_mcPartPoolIdx == m_mcPartPool.end())
    {
        // Add one more McParticle just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_mcPartPoolIdx = m_mcPartPool.insert(m_mcPartPoolIdx, McParticle());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_mcPartPool.size()), (int)MCMAXPOOLGROW);

        // For good measure expand the pool by our starting poolsize
        m_mcPartPool.insert(m_mcPartPool.end(), newSize, McParticle());
    }

    // Get the pointer to an available McParticle
    mcPart = &*m_mcPartPoolIdx++;

    return mcPart;
}

McPositionHit* McObjectManager::getNewMcPositionHit()
{
    McPositionHit* posHit = 0;

    // If we have exceeded our pre-allocated list of McParticles then expand
    if (m_mcPosHitPoolIdx == m_mcPosHitPool.end())
    {
        // Add one more McPositionHit just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_mcPosHitPoolIdx = m_mcPosHitPool.insert(m_mcPosHitPoolIdx, McPositionHit());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_mcPosHitPool.size()), (int)MCHITMAXPOOLGROW);

        // For good measure expand the pool by our starting poolsize
        m_mcPosHitPool.insert(m_mcPosHitPool.end(), newSize, McPositionHit());
    }

    // Get the pointer to an available McParticle
    posHit = &*m_mcPosHitPoolIdx++;

    return posHit;
}

McIntegratingHit* McObjectManager::getNewMcIntegratingHit()
{
    McIntegratingHit* intHit = 0;

    // If we have exceeded our pre-allocated list of McIntegratingHits then expand
    if (m_mcIntHitPoolIdx == m_mcIntHitPool.end())
    {
        // Add one more McIntegratingHit just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_mcIntHitPoolIdx = m_mcIntHitPool.insert(m_mcIntHitPoolIdx, McIntegratingHit());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_mcIntHitPool.size()), (int)MCHITMAXPOOLGROW);

        // For good measure expand the pool by our starting poolsize
        m_mcIntHitPool.insert(m_mcIntHitPool.end(), newSize, McIntegratingHit());
    }

    // Get the pointer to an available McParticle
    intHit = &*m_mcIntHitPoolIdx++;

    return intHit;
}

McTrajectory* McObjectManager::getNewMcTrajectory()
{
    McTrajectory* traj = 0;

    // If we have exceeded our pre-allocated list of McIntegratingHits then expand
    if (m_mcTrajectoryPoolIdx == m_mcTrajectoryPool.end())
    {
        // Add one more McTrajectory just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_mcTrajectoryPoolIdx = m_mcTrajectoryPool.insert(m_mcTrajectoryPoolIdx, McTrajectory());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_mcTrajectoryPool.size()), (int)MCMAXPOOLGROW);

        // For good measure expand the pool by our starting poolsize
        m_mcTrajectoryPool.insert(m_mcTrajectoryPool.end(), newSize, McTrajectory());
    }

    // Get the pointer to an available McTrajectory
    traj = &*m_mcTrajectoryPoolIdx++;

    return traj;
}

McTrajectoryPoint* McObjectManager::getNewMcTrajectoryPoint()
{
    McTrajectoryPoint* point = 0;

    // If we have exceeded our pre-allocated list of McIntegratingHits then expand
    if (m_mcTrajPointPoolIdx == m_mcTrajPointPool.end())
    {
        // Add one more McTrajectory just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_mcTrajPointPoolIdx = m_mcTrajPointPool.insert(m_mcTrajPointPoolIdx, McTrajectoryPoint());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_mcTrajPointPool.size()), (int)MCPOINTMAXPOOLGROW);

        // For good measure expand the pool by our starting poolsize
        m_mcTrajPointPool.insert(m_mcTrajPointPool.end(), newSize, McTrajectoryPoint());
    }

    // Get the pointer to an available McTrajectory
    point = &*m_mcTrajPointPoolIdx++;

    return point;
}

void McObjectManager::Delete()
{
    // Keep all of the pools down to a reasonable size...
    // Start with McParticle pool
    if (m_mcPartPool.size() > 5*MCPOOLSIZE)
    {
        m_mcPartPool.resize(5*MCPOOLSIZE);
    }

    // Now the with McPositionHit pool
    if (m_mcPosHitPool.size() > 5*MCHITPOOLSIZE)
    {
        m_mcPosHitPool.resize(5*MCHITPOOLSIZE);
    }

    // Now the with McIntegratingHit pool
    if (m_mcIntHitPool.size() > 5*MCHITPOOLSIZE)
    {
        m_mcIntHitPool.resize(5*MCHITPOOLSIZE);
    }

    // Now the with McTrajectory pool
    if (m_mcTrajectoryPool.size() > 5*MCPOOLSIZE)
    {
        m_mcTrajectoryPool.resize(5*MCPOOLSIZE);
    }

    // Now the with McTrajectoryPoint pool
    if (m_mcTrajPointPool.size() > 5*MCPOINTPOOLSIZE)
    {
        m_mcTrajPointPool.resize(5*MCPOINTPOOLSIZE);
    }

    // Ok, reset iterators to first element of our pools
    m_mcPartPoolIdx       = m_mcPartPool.begin();
    m_mcPosHitPoolIdx     = m_mcPosHitPool.begin();
    m_mcIntHitPoolIdx     = m_mcIntHitPool.begin();
    m_mcTrajectoryPoolIdx = m_mcTrajectoryPool.begin();
    m_mcTrajPointPoolIdx  = m_mcTrajPointPool.begin();

    return;
}



