#ifndef TKRRECON_H
#define TKRRECON_H

#include "TObject.h"
#include "TObjArray.h"
#include "TkrGamma.h"

// forward declarations
class TkrTrack;
class TkrSiCluster;

// ===================================
//          TrkRecon
// ===================================
//! Root object that contains all data associated with TKR Reconstruction.
/*! This is the top level object containing all data
associated with the Tracker Reconstruction.
Currently, the full reconstruction is represented by a
list of tracks (TkrTrack), a list of SiClusters (TkrSiClusters), 
and a pointer to a reconstructed gamma (TkrGamma).
*/

class TkrRecon : public TObject
{


private:
    //! list of tracks
    TObjArray *m_tracks; 

    //! Reconstructed gamma for this event
    TkrGamma m_gamma;
      
    //! list of TkrSiClusters
    /*! All TkrSiClusters found for this event
    whether or not they were associated with a track
    should there be a flag or a pointer to the TkrTrack
    to denote that a TkrSiCluster was used in a track?
    */
    TObjArray *m_SiClusters;

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

    //! provide access to the list of tracks
    /*! Clients must cast the objects to TkrTrack
    */
    TObjArray *getTracks(){ return m_tracks; };

    //! provide access to the list of clusters
    /*! Clients must cast the objects to TkrSiCluster
    */
    TObjArray *getSiClusters() { return m_SiClusters; };

    //! provide access to the gamma
    TkrGamma *getGamma() { return &m_gamma; };   
    //inline void setGamma(TkrGamma *gam) { m_gamma = gam; };

    //! access a specified TkrSiCluster from the TObjArray of TkrSiClusters
    TkrSiCluster* getSiCluster(int i) {return (TkrSiCluster*) m_SiClusters->At(i);};

    ClassDef(TkrRecon,2)


};

#endif
