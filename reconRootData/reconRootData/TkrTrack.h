#ifndef ROOT_TKRTRACK_H
#define ROOT_TKRTRACK_H

#include "TObject.h"
#include "TkrHitPlane.h"
/** 
 * @class TkrTrack
 * @brief Root object representing a reconstructed track
 *
 * Includes base summary information about the track
 * and an object array of information for each hit plane
 * This class should duplicate Event::TkrFitTrack TDS class
 * 
 * @author Tracy Usher
 * 
 * $Header$
 */

class TkrTrack : public TObject
{
public:
    //! Define the two interesting ends of the track
    enum TKREND{FIRSTHIT,LASTHIT,NOHIT};

    TkrTrack();

    virtual ~TkrTrack();

    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    //! Define two initialization methods
    void     initializeInfo(UInt_t id, UInt_t xgaps, UInt_t ygaps, 
        UInt_t x1st, UInt_t y1st);
    void     initializeQual(Double_t chiSq, Double_t ChiSqSmooth, 
        Double_t rms, Double_t quality, Double_t e, Double_t ms);

    //! Allow hits to be added to the list
    void     addHit(const TkrHitPlane& hit) {m_hits.push_back(hit);}

    //! Implement methods parallel to TkrRecInfo here
    Double_t   getQuality() const            {return m_q;}
    Double_t   getEnergy(TKREND end) const;
    UInt_t     getLayer(TKREND end) const;
    UInt_t     getTower(TKREND end) const;
    const TkrParams* getTrackPar(TKREND end) const;
    Double_t   getTrackParZ(TKREND end) const;
    const TkrCovMat* getTrackCov(TKREND end) const;

    //! Group together the access methods here
    Int_t    getId() const                   {return m_id;            }
    Double_t  getChiSq() const                {return m_chiSq;         }
    Double_t  getChiSqSmooth() const          {return m_chiSqSmooth;   }
    Double_t  getRmsResid() const             {return m_rmsResid;      }
    Double_t  getKalEnergy() const            {return m_kalEnergy;     }
    Double_t  getKalThetaMS() const           {return m_kalThetaMS;    }
    UInt_t    getXgaps() const                {return m_xGaps;         }
    UInt_t    getYgaps() const                {return m_yGaps;         }
    UInt_t    getXistGaps() const             {return m_XistGaps;      }
    UInt_t    getYistGaps() const             {return m_YistGaps;      }

    //! provide access to info at one end of the track or the other
    const TkrHitPlane* getEndHit(TKREND end) const;

    //! Provide info to any hit on the track
    UInt_t             getNumHits() const         {return m_hits.size();}
    const TkrHitPlane* getHitPlane(int idx) const {return &m_hits[idx];}

private:
    //! Summary track quality information
    Double_t   m_chiSq;
    Double_t   m_chiSqSmooth;
    Double_t   m_rmsResid;
    Double_t   m_q;

    //! Information from the Kalman Fitter
    Double_t   m_kalEnergy;
    Double_t   m_kalThetaMS;

    //! Track ID 
    Int_t     m_id;

    //! Hit Gap information
    UInt_t     m_xGaps;
    UInt_t     m_yGaps;
    UInt_t     m_XistGaps;
    UInt_t     m_YistGaps;

    //! Object array of hit planes
    //TkrHitPlaneVector m_hits;
    std::vector<TkrHitPlane> m_hits;
    ClassDef(TkrTrack,1)

};
#endif
