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
#define VECPOINTPOOLSIZE         500
#define VECPOINTMAXPOOLGROW      1000
#define VECPOINTSLINKPOOLSIZE    1000
#define VECPOINTSLINKMAXPOOLGROW 1000
#define VECNODEPOOLSIZE          500
#define VECNODEMAXPOOLGROW       500
#define TREEPOOLSIZE             20
#define TREEMAXPOOLGROW          10
#define FILTERPARAMSPOOLSIZE     20
#define FILTERPARAMSMAXPOOLGROW  10


ReconObjectManager::ReconObjectManager() : m_tkrClusterPool(CLUSTERPOOLSIZE),
                                           m_tkrTrackPool(TRACKPOOLSIZE),
                                           m_tkrTrackHitPool(TRACKHITPOOLSIZE),
                                           m_tkrTrackParamsPool(TRACKPARAMSPOOLSIZE),
                                           m_tkrVertexPool(VERTEXPOOLSIZE),
                                           m_tkrTruncationPool(TRUNCATIONPOOLSIZE),
                                           m_tkrVecPointPool(VECPOINTPOOLSIZE),
                                           m_tkrVecPointsLinkPool(VECPOINTSLINKPOOLSIZE),
                                           m_tkrVecNodePool(VECNODEPOOLSIZE),
                                           m_tkrVecNodeCompressedPool(VECNODEPOOLSIZE),
                                           m_tkrTreePool(TREEPOOLSIZE),
                                           m_tkrTreeCompressedPool(TREEPOOLSIZE),
                                           m_tkrFilterParamsPool(FILTERPARAMSPOOLSIZE)
{
    // Set pool iterators to the first item in each pool
    m_tkrClusterPoolIdx           = m_tkrClusterPool.begin();
    m_tkrTrackPoolIdx             = m_tkrTrackPool.begin();
    m_tkrTrackHitPoolIdx          = m_tkrTrackHitPool.begin();
    m_tkrTrackParamsPoolIdx       = m_tkrTrackParamsPool.begin();
    m_tkrVertexPoolIdx            = m_tkrVertexPool.begin();
    m_tkrTruncationPoolIdx        = m_tkrTruncationPool.begin();
    m_tkrVecPointPoolIdx          = m_tkrVecPointPool.begin();
    m_tkrVecPointsLinkPoolIdx     = m_tkrVecPointsLinkPool.begin();
    m_tkrVecNodePoolIdx           = m_tkrVecNodePool.begin();
    m_tkrVecNodeCompressedPoolIdx = m_tkrVecNodeCompressedPool.begin();
    m_tkrTreePoolIdx              = m_tkrTreePool.begin();
    m_tkrTreeCompressedPoolIdx    = m_tkrTreeCompressedPool.begin();
    m_tkrFilterParamsPoolIdx      = m_tkrFilterParamsPool.begin();
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

TkrVecPoint* ReconObjectManager::getNewTkrVecPoint()
{
    TkrVecPoint* vecPoint = 0;

    // If we have exceeded our pre-allocated list of TkrTrackHits then expand
    if (m_tkrVecPointPoolIdx == m_tkrVecPointPool.end())
    {
        // Add one more TkrVecPoints just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrVecPointPoolIdx = m_tkrVecPointPool.insert(m_tkrVecPointPoolIdx, TkrVecPoint());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrVecPointPool.size()), (int)VECPOINTMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrVecPointPool.insert(m_tkrVecPointPool.end(), newSize, TkrVecPoint());
    }

    // Get the pointer to an available AcdDigi
    vecPoint = &*m_tkrVecPointPoolIdx++;

    return vecPoint;
}

TkrVecPointsLink* ReconObjectManager::getNewTkrVecPointsLink()
{
    TkrVecPointsLink* vecPointsLink = 0;

    // If we have exceeded our pre-allocated list of TkrTrackHits then expand
    if (m_tkrVecPointsLinkPoolIdx == m_tkrVecPointsLinkPool.end())
    {
        // Add one more TkrVecPointsLinks just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrVecPointsLinkPoolIdx = m_tkrVecPointsLinkPool.insert(m_tkrVecPointsLinkPoolIdx, TkrVecPointsLink());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrVecPointsLinkPool.size()), (int)VECPOINTSLINKMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrVecPointsLinkPool.insert(m_tkrVecPointsLinkPool.end(), newSize, TkrVecPointsLink());
    }

    // Get the pointer to an available AcdDigi
    vecPointsLink = &*m_tkrVecPointsLinkPoolIdx++;

    return vecPointsLink;
}

TkrVecNode* ReconObjectManager::getNewTkrVecNode()
{
    TkrVecNode* vecNode = 0;

    // If we have exceeded our pre-allocated list of TkrTrackHits then expand
    if (m_tkrVecNodePoolIdx == m_tkrVecNodePool.end())
    {
        // Add one more TkrVecNodes just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrVecNodePoolIdx = m_tkrVecNodePool.insert(m_tkrVecNodePoolIdx, TkrVecNode());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrVecNodePool.size()), (int)VECNODEMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrVecNodePool.insert(m_tkrVecNodePool.end(), newSize, TkrVecNode());
    }

    // Get the pointer to an available AcdDigi
    vecNode = &*m_tkrVecNodePoolIdx++;

    return vecNode;
}

TkrVecNodeCompressed* ReconObjectManager::getNewTkrVecNodeCompressed()
{
    TkrVecNodeCompressed* vecNode = 0;

    // If we have exceeded our pre-allocated list of TkrTrackHits then expand
    if (m_tkrVecNodeCompressedPoolIdx == m_tkrVecNodeCompressedPool.end())
    {
        // Add one more TkrVecNodeCompresseds just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrVecNodeCompressedPoolIdx = m_tkrVecNodeCompressedPool.insert(m_tkrVecNodeCompressedPoolIdx, TkrVecNodeCompressed());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrVecNodeCompressedPool.size()), (int)VECNODEMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrVecNodeCompressedPool.insert(m_tkrVecNodeCompressedPool.end(), newSize, TkrVecNodeCompressed());
    }

    // Get the pointer to an available AcdDigi
    vecNode = &*m_tkrVecNodeCompressedPoolIdx++;

    return vecNode;
}

TkrTreeCompressed* ReconObjectManager::getNewTkrTreeCompressed()
{
    TkrTreeCompressed* tree = 0;

    // If we have exceeded our pre-allocated list of TkrTrackHits then expand
    if (m_tkrTreeCompressedPoolIdx == m_tkrTreeCompressedPool.end())
    {
        // Add one more TkrTrees just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrTreeCompressedPoolIdx = m_tkrTreeCompressedPool.insert(m_tkrTreeCompressedPoolIdx, TkrTreeCompressed());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrTreeCompressedPool.size()), (int)TREEMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrTreeCompressedPool.insert(m_tkrTreeCompressedPool.end(), newSize, TkrTreeCompressed());
    }

    // Get the pointer to an available AcdDigi
    tree = &*m_tkrTreeCompressedPoolIdx++;

    return tree;
}

TkrTree* ReconObjectManager::getNewTkrTree()
{
    TkrTree* tree = 0;

    // If we have exceeded our pre-allocated list of TkrTrackHits then expand
    if (m_tkrTreePoolIdx == m_tkrTreePool.end())
    {
        // Add one more TkrTrees just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrTreePoolIdx = m_tkrTreePool.insert(m_tkrTreePoolIdx, TkrTree());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrTreePool.size()), (int)TREEMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrTreePool.insert(m_tkrTreePool.end(), newSize, TkrTree());
    }

    // Get the pointer to an available AcdDigi
    tree = &*m_tkrTreePoolIdx++;

    return tree;
}

TkrFilterParams* ReconObjectManager::getNewTkrFilterParams()
{
    TkrFilterParams* filterParams = 0;

    // If we have exceeded our pre-allocated list of TkrTrackHits then expand
    if (m_tkrFilterParamsPoolIdx == m_tkrFilterParamsPool.end())
    {
        // Add one more TkrFilterParamss just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrFilterParamsPoolIdx = m_tkrFilterParamsPool.insert(m_tkrFilterParamsPoolIdx, TkrFilterParams());

        // Expand the pool by some reasonable amount
        int newSize = std::min((int)(2 * m_tkrFilterParamsPool.size()), (int)FILTERPARAMSMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrFilterParamsPool.insert(m_tkrFilterParamsPool.end(), newSize, TkrFilterParams());
    }

    // Get the pointer to an available AcdDigi
    filterParams = &*m_tkrFilterParamsPoolIdx++;

    return filterParams;
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
        m_tkrVecPointPool.clear();
        m_tkrVecPointsLinkPool.clear();
        m_tkrVecNodePool.clear();
        m_tkrVecNodeCompressedPool.clear();
        m_tkrTreePool.clear();
        m_tkrTreeCompressedPool.clear();
        m_tkrFilterParamsPool.clear();

        m_tkrClusterPoolIdx           = m_tkrClusterPool.begin();
        m_tkrTrackPoolIdx             = m_tkrTrackPool.begin();
        m_tkrTrackHitPoolIdx          = m_tkrTrackHitPool.begin();
        m_tkrTrackParamsPoolIdx       = m_tkrTrackParamsPool.begin();
        m_tkrVertexPoolIdx            = m_tkrVertexPool.begin();
        m_tkrTruncationPoolIdx        = m_tkrTruncationPool.begin();
        m_tkrVecPointPoolIdx          = m_tkrVecPointPool.begin();
        m_tkrVecPointsLinkPoolIdx     = m_tkrVecPointsLinkPool.begin();
        m_tkrVecNodePoolIdx           = m_tkrVecNodePool.begin();
        m_tkrVecNodeCompressedPoolIdx = m_tkrVecNodeCompressedPool.begin();
        m_tkrTreePoolIdx              = m_tkrTreePool.begin();
        m_tkrTreeCompressedPoolIdx    = m_tkrTreeCompressedPool.begin();
        m_tkrFilterParamsPoolIdx      = m_tkrFilterParamsPool.begin();

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

    // Now the TkrVecPoint pool
    if (m_tkrVecPointPool.size() > 5*VECPOINTPOOLSIZE)
    {
        m_tkrVecPointPool.resize(5*VECPOINTPOOLSIZE);
    }

    // Now the TkrVecPointsLink pool
    if (m_tkrVecPointsLinkPool.size() > 5*VECPOINTSLINKPOOLSIZE)
    {
        m_tkrVecPointsLinkPool.resize(5*VECPOINTSLINKPOOLSIZE);
    }

    // Now the TkrVecNode pool
    if (m_tkrVecNodePool.size() > 5*VECNODEPOOLSIZE)
    {
        m_tkrVecNodePool.resize(5*VECNODEPOOLSIZE);
    }

    // Now the TkrVecNode pool
    if (m_tkrVecNodeCompressedPool.size() > 5*VECNODEPOOLSIZE)
    {
        m_tkrVecNodeCompressedPool.resize(5*VECNODEPOOLSIZE);
    }

    // Now the TkrTree pool
    if (m_tkrTreePool.size() > 5*TREEPOOLSIZE)
    {
        m_tkrTreePool.resize(5*TREEPOOLSIZE);
    }

    // Now the TkrTree pool
    if (m_tkrTreeCompressedPool.size() > 5*TREEPOOLSIZE)
    {
        m_tkrTreeCompressedPool.resize(5*TREEPOOLSIZE);
    }

    // Now the TkrFilterParams pool
    if (m_tkrFilterParamsPool.size() > 5*FILTERPARAMSPOOLSIZE)
    {
        m_tkrFilterParamsPool.resize(5*FILTERPARAMSPOOLSIZE);
    }

    // Ok, reset iterators to first element of our pools
    m_tkrClusterPoolIdx           = m_tkrClusterPool.begin();
    m_tkrTrackPoolIdx             = m_tkrTrackPool.begin();
    m_tkrTrackHitPoolIdx          = m_tkrTrackHitPool.begin();
    m_tkrTrackParamsPoolIdx       = m_tkrTrackParamsPool.begin();
    m_tkrVertexPoolIdx            = m_tkrVertexPool.begin();
    m_tkrTruncationPoolIdx        = m_tkrTruncationPool.begin();
    m_tkrVecPointPoolIdx          = m_tkrVecPointPool.begin();
    m_tkrVecPointsLinkPoolIdx     = m_tkrVecPointsLinkPool.begin();
    m_tkrVecNodePoolIdx           = m_tkrVecNodePool.begin();
    m_tkrVecNodeCompressedPoolIdx = m_tkrVecNodeCompressedPool.begin();
    m_tkrTreePoolIdx              = m_tkrTreePool.begin();
    m_tkrTreeCompressedPoolIdx    = m_tkrTreeCompressedPool.begin();
    m_tkrFilterParamsPoolIdx      = m_tkrFilterParamsPool.begin();

    return;
}



