#ifndef ROOT_TKRRECON_H
#define ROOT_TKRRECON_H

#include "TObject.h"
#include "TObjArray.h"

#include "TkrCluster.h"
#include "TkrCandTrack.h"
#include "TkrTrack.h"
#include "TkrVertex.h"
/** 
 * @class TkrRecon
 * @brief Root container for the reconstructed Tracker information 
 * 
 * This class contains:
 * -# Collection of Clusters
 * -# Collection of 
 * -# Collectoin of Reconstructed Tracks
 * -# Collection Vertices
 * 
 * @author Tracy Usher
 *
 * $Header$
 */

class TkrRecon : public TObject
{

public:
    TkrRecon();

    virtual ~TkrRecon();

    void initialize();
    
    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    /// provide access to the list of clusters
    TObjArray *getClusterCol() { return m_clusterCol; };
    void addCluster(TkrCluster *cluster) { m_clusterCol->Add(cluster); };

    TObjArray *getTrackCandCol() { return m_trackCandCol; };
    void addTrackCand(TkrCandTrack *track) { m_trackCandCol->Add(track); };

    TObjArray *getTrackCol() { return m_trackCol; };
    void addTrack(TkrTrack* track) { m_trackCol->Add(track); };

    TObjArray *getVertexCol() { return m_vertexCol; };
    void addVertex(TkrVertex* vertex) { m_vertexCol->Add(vertex); };

private:
    /// collection of TkrSiCluster
    TObjArray *m_clusterCol;
    /// collection of TkrCandTrack
    TObjArray *m_trackCandCol;
    /// collection of TkrTrack
    TObjArray *m_trackCol;
    /// collection of TkrVertex
    TObjArray *m_vertexCol;

    ClassDef(TkrRecon,4)
};

#endif
