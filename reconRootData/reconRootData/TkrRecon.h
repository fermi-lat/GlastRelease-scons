#ifndef TKRRECON_H
#define TKRRECON_H

#include "TObject.h"
#include "TObjArray.h"

#include "TkrCandTrack.h"
#include "TkrSiCluster.h"
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
    
    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    /// provide access to the list of clusters
    /// Clients must cast the objects to TkrSiCluster
    TObjArray* getSiClusterCol() {return m_siClusterCol;};
    void addSiCluster(TkrSiCluster* cluster) { m_siClusterCol->Add(cluster); };

    /// provide access to the list of pattern recognized tracks
    /// Clients must cast the objects to TkrTrackCand
    TObjArray* getTrackCandCol() {return m_trackCandCol;};
    void addTrackCand(TkrCandTrack *track) { m_trackCandCol->Add(track); };

    /// provide access to the list of tracks
    /// Clients must cast the objects to TkrTrack
    TObjArray* getTrackCol() {return m_trackCol;};
    void addTrack(TkrTrack *track) { m_trackCol->Add(track); };

    /// provide access to the list of vertices
    /// Clients must cast the objects to TkrVertex
    TObjArray* getVertexCol() {return m_vertexCol;};
    void addVertex(TkrVertex *vertex) { m_vertexCol->Add(vertex); };


private:
    //! list of TkrSiClusters
    TObjArray* m_siClusterCol;
    static TObjArray *s_siClusterStaticCol; //! Collection of TkrSiCluster objects

    //! list of TkrTrackCands
    TObjArray* m_trackCandCol;
    static TObjArray *s_trackCandStaticCol; //! Collection of TkrCandTrack objects

    //! list of TkrTracks
    TObjArray* m_trackCol;
    static TObjArray *s_trackStaticCol; //! Collection of TkrTrack objects

    //! list of TkrVertices
    TObjArray* m_vertexCol;
    static TObjArray *s_vertexStaticCol; //! Collection of TkrVertex objects


    ClassDef(TkrRecon,3)
};

#endif
