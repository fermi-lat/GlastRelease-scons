#ifndef ROOT_TKRTRACK_H
#define ROOT_TKRTRACK_H

#include "TObject.h"
#include "TObjArray.h"
#include "TkrHitPlane.h"
/** 
 * @class TkrTrack
 *
 * @brief Root object representing a reconstructed track
 *        Includes base summary information about the track
 *        and an object array of information for each hit plane
 *        This class should duplicate TkrFitTrack
 * 
 * Created 13-May-2002
 * 
 * @author Tracy Usher
 *
 */

class TkrTrack : public TObject
{
public:
    //! Define the two interesting ends of the track
    enum TKREND{FIRSTHIT,LASTHIT,NOHIT};

    TkrTrack();

    virtual ~TkrTrack() {m_hits.clear();}

    //! Define two initialization methods
    void     initializeInfo(UInt_t id, UInt_t xgaps, UInt_t ygaps, UInt_t x1st, UInt_t y1st);
    void     initializeQaul(Double_t chiSq, Double_t ChiSqSmooth, Double_t rms, Double_t quality, Double_t e, Double_t ms);

    //! Allow hits to be added to the list
    void     addHit(const TkrHitPlane& hit) {m_hits.push_back(hit);}

    //! Implement methods parallel to TkrRecInfo here
    const Double_t   getQuality()             {return m_Q;                 }
    const Double_t   getEnergy(TKREND end)    {getEndHit(end).getEnePlane(); }
    const UInt_t     getLayer(TKREND end)     {getEndHit(end).getIdPlane();}
    const UInt_t     getTower(TKREND end)     {getEndHit(end).getIdTower();}
    const TkrParams& getTrackPar(TKREND end)  {getEndHit(end).getHitSmooth().getTkrParams();}
    const Double_t   getTrackParZ(TKREND end) {getEndHit(end).getZplane();}
    const TkrCovMat& getTrackCov(TKREND end)  {getEndHit(end).getHitSmooth().getTkrCovMat();}

    //! Group together the access methods here
    UInt_t    getId()                         {return m_id;            }
    Double_t  getChiSq()                      {return m_ChiSq;         }
    Double_t  getChiSqSmooth()                {return m_ChiSqSmooth;   }
    Double_t  getRmsResid()                   {return m_rmsResid;      }
    Double_t  getQ()                          {return m_Q;             }
    Double_t  getKalEnergy()                  {return m_KalEnergy;     }
    Double_t  getKalThetaMS()                 {return m_KalThetaMS;    }
    UInt_t    getXgaps()                      {return m_Xgaps;         }
    UInt_t    getYgaps()                      {return m_Ygaps;         }
    UInt_t    getXistGaps()                   {return m_XistGaps;      }
    UInt_t    getYistGaps()                   {return m_YistGaps;      }

    //! provide access to info at one end of the track or the other
    const TkrHitPlane& getEndHit(TKREND end);

    //! Provide info to any hit on the track
    UInt_t             getNumHits()           {return m_hits.size();}
    TkrHitPlaneIter    getHitIterBegin()      {return m_hits.begin();  }
    TkrHitPlaneIter    getHitIterEnd()        {return m_hits.end();    }

private:
    //! Track ID 
    UInt_t     m_id;

    //! Summary track quality information
    Double_t   m_ChiSq;
    Double_t   m_ChiSqSmooth;
    Double_t   m_rmsResid;
    Double_t   m_Q;

    //! Information from the Kalman Fitter
    Double_t   m_KalEnergy;
    Double_t   m_KalThetaMS;

    //! Hit Gap information
    UInt_t     m_Xgaps;
    UInt_t     m_Ygaps;
    UInt_t     m_XistGaps;
    UInt_t     m_YistGaps;

    //! Object array of hit planes
    TkrHitPlaneVector m_hits;

    ClassDef(TkrTrack,1)

};
#endif
