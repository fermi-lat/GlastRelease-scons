// File and Version Information:
// $Header$
//
// Description: this utility singleton is used in various classes of G4Generator
// to register new McParticle objects, retrive the actual McParticle (i.e. the
// one that is actually creating hits in the sensitive detectors) and to finally
// save the McParticle hierarchy in the TDS
//      
// Author(s):
//      Iama Crazy

#include "ReconObjectManager.h"

// This is the singleton static pointer
ReconObjectManager* ReconObjectManager::m_pointer = 0;

// Define the default pool sizes for our objects
#define CLUSTERPOOLSIZE          100
#define CLUSTERMAXPOOLGROW       1000
#define TRACKPOOLSIZE            50
#define TRACKMAXPOOLGROW         100
#define TRACKHITPOOLSIZE         500
#define TRACKHITMAXPOOLGROW      1000
#define TRACKPARAMSPOOLSIZE      500
#define TRACKPARAMSMAXPOOLGROW   1000
#define VERTEXPOOLSIZE           20
#define VERTEXMAXPOOLGROW        50
#define TRUNCATIONPOOLSIZE       50
#define TRUNCATIONMAXPOOLGROW    100

ReconObjectManager::ReconObjectManager() : m_tkrClusterPool(CLUSTERPOOLSIZE),
                                           m_tkrTrackPool(TRACKPOOLSIZE),
                                           m_tkrTrackHitPool(TRACKHITPOOLSIZE),
                                           m_tkrTrackParamsPool(TRACKPARAMSPOOLSIZE),
                                           m_tkrVertexPool(VERTEXPOOLSIZE),
                                           m_tkrTruncationPool(TRUNCATIONPOOLSIZE)
{
    // Set pool iterators to the first item in each pool
    m_tkrClusterPoolIdx     = m_tkrClusterPool.begin();
    m_tkrTrackPoolIdx       = m_tkrTrackPool.begin();
    m_tkrTrackHitPoolIdx    = m_tkrTrackHitPool.begin();
    m_tkrTrackParamsPoolIdx = m_tkrTrackParamsPool.begin();
    m_tkrVertexPoolIdx      = m_tkrVertexPool.begin();
    m_tkrTruncationPoolIdx  = m_tkrTruncationPool.begin();
}

ReconObjectManager* ReconObjectManager::getPointer()
{
  // Purpose and Method: standard singleton method to retrive the unique pointer
  if(m_pointer == 0) m_pointer = new ReconObjectManager();
  return m_pointer;
}

TkrCluster* ReconObjectManager::getNewTkrCluster()
{
    TkrCluster* cluster = 0;

    // If we have exceeded our pre-allocated list of TkrClusters then expand
    if (m_tkrClusterPoolIdx == m_tkrClusterPool.end())
    {
        // Add one more TkrCluster just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrClusterPoolIdx = m_tkrClusterPool.insert(m_tkrClusterPoolIdx, TkrCluster());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrClusterPool.size()), (int)CLUSTERMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrClusterPool.insert(m_tkrClusterPool.end(), newSize, TkrCluster());
    }

    // Get the pointer to an available AcdDigi
    cluster = &*m_tkrClusterPoolIdx++;

    return cluster;
}

TkrTrack* ReconObjectManager::getNewTkrTrack()
{
    TkrTrack* track = 0;

    // If we have exceeded our pre-allocated list of TkrTracks then expand
    if (m_tkrTrackPoolIdx == m_tkrTrackPool.end())
    {
        // Add one more TkrTrack just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrTrackPoolIdx = m_tkrTrackPool.insert(m_tkrTrackPoolIdx, TkrTrack());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrTrackPool.size()), (int)TRACKMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrTrackPool.insert(m_tkrTrackPool.end(), newSize, TkrTrack());
    }

    // Get the pointer to an available AcdDigi
    track = &*m_tkrTrackPoolIdx++;

    return track;
}

TkrTrackHit* ReconObjectManager::getNewTkrTrackHit()
{
    TkrTrackHit* trackHit = 0;

    // If we have exceeded our pre-allocated list of TkrTrackHits then expand
    if (m_tkrTrackHitPoolIdx == m_tkrTrackHitPool.end())
    {
        // Add one more TkrTrackHits just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrTrackHitPoolIdx = m_tkrTrackHitPool.insert(m_tkrTrackHitPoolIdx, TkrTrackHit());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrTrackHitPool.size()), (int)TRACKHITMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrTrackHitPool.insert(m_tkrTrackHitPool.end(), newSize, TkrTrackHit());
    }

    // Get the pointer to an available AcdDigi
    trackHit = &*m_tkrTrackHitPoolIdx++;

    return trackHit;
}

TkrTrackParams* ReconObjectManager::getNewTkrTrackParams()
{
    TkrTrackParams* trackHit = 0;

    // If we have exceeded our pre-allocated list of TkrTrackParamss then expand
    if (m_tkrTrackParamsPoolIdx == m_tkrTrackParamsPool.end())
    {
        // Add one more TkrTrackParamss just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrTrackParamsPoolIdx = m_tkrTrackParamsPool.insert(m_tkrTrackParamsPoolIdx, TkrTrackParams());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrTrackParamsPool.size()), (int)TRACKPARAMSMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrTrackParamsPool.insert(m_tkrTrackParamsPool.end(), newSize, TkrTrackParams());
    }

    // Get the pointer to an available AcdDigi
    trackHit = &*m_tkrTrackParamsPoolIdx++;

    return trackHit;
}

TkrVertex* ReconObjectManager::getNewTkrVertex()
{
    TkrVertex* vertex = 0;

    // If we have exceeded our pre-allocated list of TkrVertexs then expand
    if (m_tkrVertexPoolIdx == m_tkrVertexPool.end())
    {
        // Add one more TkrVertexs just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrVertexPoolIdx = m_tkrVertexPool.insert(m_tkrVertexPoolIdx, TkrVertex());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrVertexPool.size()), (int) VERTEXMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrVertexPool.insert(m_tkrVertexPool.end(), newSize, TkrVertex());
    }

    // Get the pointer to an available AcdDigi
    vertex = &*m_tkrVertexPoolIdx++;

    return vertex;
}

TkrTruncationData* ReconObjectManager::getNewTkrTruncationData()
{
    TkrTruncationData* truncation = 0;

    // If we have exceeded our pre-allocated list of TkrVertexs then expand
    if (m_tkrTruncationPoolIdx == m_tkrTruncationPool.end())
    {
        // Add one more TkrVertexs just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrTruncationPoolIdx = m_tkrTruncationPool.insert(m_tkrTruncationPoolIdx, TkrTruncationData());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrTruncationPool.size()), (int) TRUNCATIONMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrTruncationPool.insert(m_tkrTruncationPool.end(), newSize, TkrTruncationData());
    }

    // Get the pointer to an available AcdDigi
    truncation = &*m_tkrTruncationPoolIdx++;

    return truncation;
}

void ReconObjectManager::Delete(const char* opt)
{
    if (opt[0] == 'A') {

        m_tkrClusterPool.clear();
        m_tkrTrackPool.clear();
        m_tkrTrackHitPool.clear();
        m_tkrTrackParamsPool.clear();
        m_tkrVertexPool.clear();
        m_tkrTruncationPool.clear();

        m_tkrClusterPoolIdx     = m_tkrClusterPool.begin();
        m_tkrTrackPoolIdx       = m_tkrTrackPool.begin();
        m_tkrTrackHitPoolIdx    = m_tkrTrackHitPool.begin();
        m_tkrTrackParamsPoolIdx = m_tkrTrackParamsPool.begin();
        m_tkrVertexPoolIdx      = m_tkrVertexPool.begin();
        m_tkrTruncationPoolIdx  = m_tkrTruncationPool.begin();
        return;
    }

    // Keep all of the pools down to a reasonable size...
    // Start with TkrCluster pool
    if (m_tkrClusterPool.size() > 5*CLUSTERPOOLSIZE)
    {
        m_tkrClusterPool.resize(5*CLUSTERPOOLSIZE);
    }

    // Now the TkrTrack pool
    if (m_tkrTrackPool.size() > 5*TRACKPOOLSIZE)
    {
        m_tkrTrackPool.resize(5*TRACKPOOLSIZE);
    }

    // Now the TkrTrackHit pool
    if (m_tkrTrackHitPool.size() > 5*TRACKHITPOOLSIZE)
    {
        m_tkrTrackHitPool.resize(5*TRACKHITPOOLSIZE);
    }

    // Now the TkrTrackParams pool
    if (m_tkrTrackParamsPool.size() > 5*TRACKPARAMSPOOLSIZE)
    {
        m_tkrTrackParamsPool.resize(5*TRACKPARAMSPOOLSIZE);
    }

    // Now the TkrVertex pool
    if (m_tkrVertexPool.size() > 5*VERTEXPOOLSIZE)
    {
        m_tkrVertexPool.resize(5*VERTEXPOOLSIZE);
    }

    // Now the TkrTruncation pool
    if (m_tkrTruncationPool.size() > 5*TRUNCATIONPOOLSIZE)
    {
        m_tkrTruncationPool.resize(5*TRUNCATIONPOOLSIZE);
    }

    // Ok, reset iterators to first element of our pools
    m_tkrClusterPoolIdx     = m_tkrClusterPool.begin();
    m_tkrTrackPoolIdx       = m_tkrTrackPool.begin();
    m_tkrTrackHitPoolIdx    = m_tkrTrackHitPool.begin();
    m_tkrTrackParamsPoolIdx = m_tkrTrackParamsPool.begin();
    m_tkrVertexPoolIdx      = m_tkrVertexPool.begin();
    m_tkrTruncationPoolIdx  = m_tkrTruncationPool.begin();
    return;
}



