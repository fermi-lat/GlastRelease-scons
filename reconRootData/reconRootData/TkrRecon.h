#ifndef TKRRECON_H
#define TKRRECON_H

#include "TObject.h"
#include "TObjArray.h"
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
    // constructors
    //! default constructor
    /*! clients must call Create() as necessary
    */
    TkrRecon();

    //! destructor
    virtual ~TkrRecon();
    
    //! release all dynamically allocated memory
    void Clean(); 

    //! Allocate memory as necessary for this class
    void Create();

    //! provide access to the list of clusters
    /*! Clients must cast the objects to TkrSiCluster
    */
    TObjArray*    getSiClusters()     {return m_SiClusters;};

    //! provide access to the list of pattern recognized tracks
    /*! Clients must cast the objects to TkrTrackCand
    */
    TObjArray*    getTrackCands()     {return m_TrackCands;};

    //! provide access to the list of tracks
    /*! Clients must cast the objects to TkrTrack
    */
    TObjArray*    getTracks()         {return m_Tracks;};

    //! provide access to the list of vertices
    /*! Clients must cast the objects to TkrVertex
    */
    TObjArray*    getVertices()       {return m_Vertices;};

    ClassDef(TkrRecon,2)

private:
    //! list of TkrSiClusters
    TObjArray* m_SiClusters;
    
    //! list of TkrTrackCands
    TObjArray* m_TrackCands;

    //! list of TkrTracks
    TObjArray* m_Tracks;

    //! list of TkrVertices
    TObjArray* m_Vertices;
};

#endif
