#ifndef TKRTRACK_H
#define TKRTRACK_H

#include "TObject.h"
#include "TObjArray.h"

// ===================================
//          TkrTrack
// ===================================
//  
//! Root object representing a reconstructed particle track.
/*! Represents a reconstructed particle track. It 
could be a charged particle track or a photon that is 
reconstructed due to it's conversion into an 
electron-positron pair. 

A track contains two lists:

1.) A list of TkrLocator objects which contain the position and angle of the 
    track at different interesting positions  (e.g.: initial point, exit point 
    of the tracker).

2.) A list of TkrHit objects which point to the TkrSiCluster objects that are associated with this track.

*/

class TkrTrack : public TObject
{
public:

    //! our set of particles types
    typedef enum {
        NONE = -1,
        GAMMA = 0,
        ELECTRON,
        POSITRON,
        HADRON,
        OTHER
    } ParticleType;

private:
    //! list of hits associated with this track
    TObjArray *m_hits;          

    //! list of TkrLocators
    TObjArray *m_locator;       
    
    //! track id - probably the position in a TObjArray?
    UShort_t m_id;        
    
    //! type of the track (gamma, e+, e-, hadron ...)
    ParticleType m_type;        

    //! id of the mother particle (-1 if none)
    ParticleType m_motherId;    

    //! residual  
    Float_t m_residual;         

    //! chi**2 of the fit
    Float_t m_chi2;             

    //! quality measure of the track
    Float_t m_quality;          

    //! layer where the particle starts
    UInt_t m_firstLayer;         

    //! energy input from the cal
    Float_t m_energyInput;    
    
    //! energy determined using the tracker + cal
    Float_t m_energyDetermined;    

public:

    //! constructors
    TkrTrack();
    TkrTrack(UShort_t id);

    //! destructor
    virtual ~TkrTrack();

    //! frees all dynamically allocated memory
    /*!  Clients should call Clean() before calling Create() */
    void Clean(); 

    //! allocates memory - if necessary
    /*! if the objects are already allocated, 
    no further memory is allocated */
    void Create();

    //! provide access to the list of hits
    /*! allows the caller to modify the TObjArray */
    TObjArray *getHits() { return m_hits; };

    //! provide access to the list of locators
    /*! allow the caller to modify the TObjArray */
    TObjArray *getLocator() { return m_locator; };

    //! access the ID # of this track
    UShort_t getId() { return m_id; };
    void setId(UShort_t id) { m_id = id; };

    //! access the track type
    ParticleType getTrackType() { return m_type; };
    void setTrackType(ParticleType type) { m_type = type; };

    //! access the mother id
    ParticleType getMother() { return m_motherId; };
    void setMother(ParticleType mom) { m_motherId = mom; };

    //! access the residual associated with this track
    Float_t getResidual() { return m_residual; };
    void setResidual(Float_t resid) { m_residual = resid; };

    //! access Chi2
    Float_t getChi2() { return m_chi2; };
    void setChi2(Float_t chi) { m_chi2 = chi; };

    //! access the quality parameter
    Float_t getQuality() { return m_quality; };
    void setQuality(Float_t q) { m_quality = q; };

    //! access the first layer in the TKR associated with this track
    UInt_t getFirstLayer() { return m_firstLayer; };
    void setFirstLayer(UInt_t first) { m_firstLayer = first; };
    
    //! access the number of TkrSiClusters associated with this track
    UInt_t getNumClusters() { return m_hits->GetEntries(); };

    //! access the input energy
    Float_t getEnergyInput() { return m_energyInput; };
    void setEnergyInput(Float_t ene) { m_energyInput = ene; };

    //! access the determined energy
    Float_t getEnergyDet() { return m_energyDetermined; };
    void setEnergyDet(Float_t ene) { m_energyDetermined = ene; };

    ClassDef(TkrTrack,1)
};
#endif
