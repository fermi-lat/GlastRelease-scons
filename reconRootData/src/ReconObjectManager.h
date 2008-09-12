#ifndef ReconObjectManager_H
#define ReconObjectManager_H

#include <list>

#include "reconRootData/TkrCluster.h"
#include "reconRootData/TkrTrack.h"
#include "reconRootData/TkrTrackHit.h"
#include "reconRootData/TkrTrackParams.h"
#include "reconRootData/TkrVertex.h"

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
    TkrCluster*     getNewTkrCluster();

    /// This method returns a "new" TkrTrack object allocated from the local pool
    TkrTrack*       getNewTkrTrack();

    /// This method returns a "new" TkrTrackHit object allocated from the local pool
    TkrTrackHit*    getNewTkrTrackHit();

    /// This method returns a "new" TkrTrackParams object allocated from the local pool
    TkrTrackParams* getNewTkrTrackParams();

    /// This method returns a "new" TkrVertex object allocated from the local pool
    TkrVertex*      getNewTkrVertex();

    /// "Delete" all objects
    void Delete(const char* opt="");

private:
    /// The constructor is private since this is a singleton
    ReconObjectManager();   

    /// The static pointer of the singleton
    static ReconObjectManager*          m_pointer;

    /// Define a "pool" for TkrCluster and an iterator for accessing them
    std::list<TkrCluster>               m_tkrClusterPool;
    std::list<TkrCluster>::iterator     m_tkrClusterPoolIdx;

    /// Define a "pool" for TkrTrack and an iterator for accessing them
    std::list<TkrTrack>                 m_tkrTrackPool;
    std::list<TkrTrack>::iterator       m_tkrTrackPoolIdx;

    /// Define a "pool" for TkrTrackHit and an iterator for accessing them
    std::list<TkrTrackHit>              m_tkrTrackHitPool;
    std::list<TkrTrackHit>::iterator    m_tkrTrackHitPoolIdx;

    /// Define a "pool" for TkrTrackParams and an iterator for accessing them
    std::list<TkrTrackParams>           m_tkrTrackParamsPool;
    std::list<TkrTrackParams>::iterator m_tkrTrackParamsPoolIdx;

    /// Define a "pool" for TkrVertex and an iterator for accessing them
    std::list<TkrVertex>                m_tkrVertexPool;
    std::list<TkrVertex>::iterator      m_tkrVertexPoolIdx;
};
#endif //ReconObjectManager_H
