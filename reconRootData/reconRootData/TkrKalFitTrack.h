#ifndef ROOT_TKRKALFITTRACK_H
#define ROOT_TKRKALFITTRACK_H

#include "TObject.h"
#include "TObjArray.h"
#include "TVector3.h"
#include "TkrHitPlane.h"
/** 
 * @class TkrKalFitTrack
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

class TkrKalFitTrack : public TObject
{
public:
    //! Define the two interesting ends of the track
    enum TKREND{FIRSTHIT,LASTHIT,NOHIT};
    enum Status{EMPTY,FOUND,CRACK};

    TkrKalFitTrack();

    virtual ~TkrKalFitTrack();

    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    //! Define four initialization methods
    void     initializeBase(Status stat, UInt_t type, 
        UInt_t id, Double_t energy0, TVector3& x0, TVector3& dir);
    void     initializeQual(Double_t chiSq, Double_t ChiSqSmooth, 
        Double_t rms, Double_t quality, Double_t e, Double_t e_err, Double_t ms);
    void     initializeGaps(UInt_t xgaps, UInt_t ygaps, UInt_t x1st, 
        UInt_t y1st);
    void     initializeKal(UInt_t nSegPoints, Double_t chisqSeg, UInt_t nxHits, 
        UInt_t nyHits, Double_t radlen);

    //! Allow hits to be added to the list
    void     addHit(const TkrHitPlane& hit) {m_hits.push_back(hit);}

    //! Implement methods parallel to Event::TkrKalFitTrackBase here
    Double_t   getQuality() const            {return m_Q;}
    Double_t   getEnergy(TKREND end) const;
    UInt_t     getLayer(TKREND end) const;
    UInt_t     getTower(TKREND end) const;
    const TkrParams* getTrackPar(TKREND end) const;
    Double_t   getTrackParZ(TKREND end) const;
    const TkrCovMat* getTrackCov(TKREND end) const;

    //! Implement methods parallel to Event::TkrKalFitTrack here
    Status     status() const                  {return m_status;          }
    UInt_t     getType() const                 {return m_type;            }
    Int_t      getId() const                   {return m_id;              }
    Double_t   getChiSquare() const            {return m_chisq;           }
    Double_t   getChiSquareSmooth() const      {return m_chisqSmooth;     }
    Double_t   getScatter() const              {return m_rmsResid;        }
    Double_t   getKalThetaMS() const           {return m_KalThetaMS;      }
    Double_t   getKalEnergy() const            {return m_KalEnergy;       }
    Double_t   getKalEnergyError() const       {return m_KalEnergyErr;    }
    UInt_t     getNumGaps() const              {return m_Xgaps+m_Ygaps;   }
    UInt_t     getNumXGaps() const             {return m_Xgaps;           }
    UInt_t     getNumYGaps() const             {return m_Ygaps;           }
    UInt_t     getNumXFirstGaps() const        {return m_XistGaps;        }
    UInt_t     getNumYFirstGaps() const        {return m_YistGaps;        }
    TVector3   getInitialPosition() const      {return m_x0;              }
    TVector3   getInitialDirection() const     {return m_dir;             }
    Double_t   getStartEnergy()const           {return m_energy0;         }
    UInt_t     getNumSegmentPoints() const     {return m_numSegmentPoints;}
    Double_t   getChiSquareSegment() const     {return m_chisqSegment;    }
    UInt_t     getNumXHits() const             {return m_nxHits;          }
    UInt_t     getNumYHits() const             {return m_nyHits;          }
    Double_t   getTkrCalRadlen() const         {return m_TkrCal_radlen;   }

    //! provide access to info at one end of the track or the other
    const TkrHitPlane* getEndHit(TKREND end) const;

    //! Provide info to any hit on the track
    UInt_t             getNumHits() const         {return m_hits.size();}
    const TkrHitPlane* getHitPlane(int idx) const {return &m_hits[idx];}

private:
    /// Status
    Status     m_status;
    UInt_t     m_type; 

    //! Track ID 
    UInt_t     m_id;

    /// The input energy, and current position and direction
    Double_t   m_energy0;
    TVector3   m_x0;
    TVector3   m_dir;

    // Track quality information
    Double_t   m_chisq;                // Standard track chi-square
    Double_t   m_chisqSmooth;          // "Smoothed" track chi-square
    Double_t   m_rmsResid;             // RMS residuals of track hits to fit
    Double_t   m_Q;                    // Track "Quality"

    // Information from the fitter
    Double_t   m_KalEnergy;
    Double_t   m_KalThetaMS;

    // Hit gap information
    UInt_t     m_Xgaps;
    UInt_t     m_Ygaps;
    UInt_t     m_XistGaps;
    UInt_t     m_YistGaps;

    /// KalTrack data
    UInt_t     m_numSegmentPoints;
    Double_t   m_chisqSegment;
    UInt_t     m_nxHits;
    UInt_t     m_nyHits;
    Double_t   m_KalEnergyErr;
    Double_t   m_TkrCal_radlen; 

    //! Object array of hit planes
    std::vector<TkrHitPlane> m_hits;

    ClassDef(TkrKalFitTrack,1)
};
#endif
