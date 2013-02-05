#ifndef ReconObjectManager_H
#define ReconObjectManager_H

#include <list>

#include "reconRootData/TkrCluster.h"
#include "reconRootData/TkrTrack.h"
#include "reconRootData/TkrTrackHit.h"
#include "reconRootData/TkrTrackParams.h"
#include "reconRootData/TkrVertex.h"
#include "reconRootData/TkrTruncationData.h"
#include "reconRootData/TkrTree.h"
#include "reconRootData/TkrFilterParams.h"

/** 
 *  @class ReconObjectManager
 *
 *  @brief A singleton manager for handling the RECON root objects
 *         This uses stl lists to manage a pool of objects for each root object in
 *         order to avoid using "new" and "delete". The pool can grow at a sensible 
 *         rate, it will automatically contract to a sensible size as well to try to
 *         prevent memory issues. 
 * 
 *         The objects are retrieved from the pool by overriding the "new" and "delete"
 *         operators in the specific classes. See, for example, TkrCluster.h in this package
 *
 *  ****   Currently only active for Tkr Recon objects (March 27, 2008)
 *
 *  @author Iama Crazy
 *
 * $Header$
 */
class ReconObjectManager {
public:
    virtual ~ReconObjectManager() {m_pointer = 0;};

    /// The static pointer retrival method of the singleton
    static ReconObjectManager* getPointer();

    /// This method returns a "new" TkrCluster object allocated from the local pool
    TkrCluster*         getNewTkrCluster();

    /// This method returns a "new" TkrTrack object allocated from the local pool
    TkrTrack*           getNewTkrTrack();

    /// This method returns a "new" TkrTrackHit object allocated from the local pool
    TkrTrackHit*        getNewTkrTrackHit();

    /// This method returns a "new" TkrTrackParams object allocated from the local pool
    TkrTrackParams*     getNewTkrTrackParams();

    /// This method returns a "new" TkrVertex object allocated from the local pool
    TkrVertex*          getNewTkrVertex();

    /// This method returns a "new" TkrVertex object allocated from the local pool
    TkrTruncationData*  getNewTkrTruncationData();

    /// This method returns a "new" TkrVecPoint object allocated from the local pool
    TkrVecPoint*        getNewTkrVecPoint();

    /// This method returns a "new" TkrVecPointsLink object allocated from the local pool
    TkrVecPointsLink*   getNewTkrVecPointsLink();

    /// This method returns a "new" TkrVecNode object allocated from the local pool
    TkrVecNode*         getNewTkrVecNode();

    /// This method returns a "new" TkrTree object allocated from the local pool
    TkrTree*            getNewTkrTree();

    /// This method returns a "new" TkrFilterParams object allocated from the local pool
    TkrFilterParams*    getNewTkrFilterParams();

    /// "Delete" all objects
    void Delete(const char* opt="");

private:
    /// The constructor is private since this is a singleton
    ReconObjectManager();   

    /// The static pointer of the singleton
    static ReconObjectManager*              m_pointer;

    /// Define a "pool" for TkrCluster and an iterator for accessing them
    std::list<TkrCluster>                   m_tkrClusterPool;
    std::list<TkrCluster>::iterator         m_tkrClusterPoolIdx;

    /// Define a "pool" for TkrTrack and an iterator for accessing them
    std::list<TkrTrack>                     m_tkrTrackPool;
    std::list<TkrTrack>::iterator           m_tkrTrackPoolIdx;

    /// Define a "pool" for TkrTrackHit and an iterator for accessing them
    std::list<TkrTrackHit>                  m_tkrTrackHitPool;
    std::list<TkrTrackHit>::iterator        m_tkrTrackHitPoolIdx;

    /// Define a "pool" for TkrTrackParams and an iterator for accessing them
    std::list<TkrTrackParams>               m_tkrTrackParamsPool;
    std::list<TkrTrackParams>::iterator     m_tkrTrackParamsPoolIdx;

    /// Define a "pool" for TkrVertex and an iterator for accessing them
    std::list<TkrVertex>                    m_tkrVertexPool;
    std::list<TkrVertex>::iterator          m_tkrVertexPoolIdx;

    /// Define a "pool" for TkrTruncationData and an iterator for accessing them
    std::list<TkrTruncationData>            m_tkrTruncationPool;
    std::list<TkrTruncationData>::iterator  m_tkrTruncationPoolIdx;

    /// Define a "pool" for TkrVecPoint objects and an iterator for accessing them
    std::list<TkrVecPoint>                  m_tkrVecPointPool;
    std::list<TkrVecPoint>::iterator        m_tkrVecPointPoolIdx;

    /// Define a "pool" for TkrVecPointsLink objects and an iterator for accessing them
    std::list<TkrVecPointsLink>             m_tkrVecPointsLinkPool;
    std::list<TkrVecPointsLink>::iterator   m_tkrVecPointsLinkPoolIdx;

    /// Define a "pool" for TkrVecNode objects and an iterator for accessing them
    std::list<TkrVecNode>                   m_tkrVecNodePool;
    std::list<TkrVecNode>::iterator         m_tkrVecNodePoolIdx;

    /// Define a "pool" for TkrTree objects and an iterator for accessing them
    std::list<TkrTree>                      m_tkrTreePool;
    std::list<TkrTree>::iterator            m_tkrTreePoolIdx;

    /// Define a "pool" for TkrFilterParams objects and an iterator for accessing them
    std::list<TkrFilterParams>              m_tkrFilterParamsPool;
    std::list<TkrFilterParams>::iterator    m_tkrFilterParamsPoolIdx;

};
#endif //ReconObjectManager_H
