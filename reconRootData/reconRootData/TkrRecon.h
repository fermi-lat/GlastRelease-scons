#ifndef ROOT_TKRRECON_H
#define ROOT_TKRRECON_H

#include "TObject.h"
#include "TObjArray.h"

#include "TkrSiCluster.h"
#include "TkrCandTrack.h"
#include "TkrTrack.h"
#include "TkrVertex.h"
/** 
 * @class TkrRecon
 *
 * @brief Root container for the reconstructed Tracker information 
 *        Reconstructed information includes:
 *           1) Silicon Clusters
 *           2) Pattern Recognized Tracks
 *           3) Reconstructed Tracks
 *           4) Vertices
 * 
 * Created 13-May-2002
 * 
 * @author Tracy Usher
 *
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
    TObjArray *getSiClusterCol() { return m_siClusterCol; };
    void addSiCluster(TkrSiCluster *cluster) { m_siClusterCol->Add(cluster); };

    TObjArray *getTrackCandCol() { return m_trackCandCol; };
    void addTrackCand(TkrCandTrack *track) { m_trackCandCol->Add(track); };

    TObjArray *getTrackCol() { return m_trackCol; };
    void addTrack(TkrTrack* track) { m_trackCol->Add(track); };

    TObjArray *getVertexCol() { return m_vertexCol; };
    void addVertex(TkrVertex* vertex) { m_vertexCol->Add(vertex); };

private:
    /// collection of TkrSiCluster
    TObjArray *m_siClusterCol;
    /// collection of TkrCandTrack
    TObjArray *m_trackCandCol;
    /// collection of TkrTrack
    TObjArray *m_trackCol;
    /// collection of TkrVertex
    TObjArray *m_vertexCol;

    ClassDef(TkrRecon,3)
};

#endif
