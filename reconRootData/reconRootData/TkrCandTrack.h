#ifndef ROOT_TKRCANDTRACK_H
#define ROOT_TKRCANDTRACK_H

#include "TObject.h"
#include "TkrCandHit.h"
#include "TVector3.h"

/** 
 * @class TkrCandTrack
 *
 * @brief Root object representing a pattern recognition candidate track
 *        Includes base summary information about the track
 *        and a vector of objects containing information for each hit plane
 *        This class should exactly duplicate TkrPatCand
 * 
 * Created 15-May-2002
 * 
 * @author Tracy Usher
 *
 */

class TkrCandTrack : public TObject
{
public:
    TkrCandTrack();

    TkrCandTrack(UInt_t id, UInt_t layer, UInt_t tower, Double_t q, Double_t e, const TVector3& pos, const TVector3& dir);

    virtual ~TkrCandTrack() {m_hits.clear();}

    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    //! Define initialization method
    void     initialize(UInt_t id, UInt_t layer, UInt_t tower, Double_t q, Double_t e, const TVector3& pos, const TVector3& dir);

    //! Allow hits to be added to the list
    void     addHit(const TkrCandHit& hit) {m_hits.push_back(hit);}

    //! Group together the access methods here
    UInt_t        getId()               const {return m_id;        }
    Double_t      getQuality()          const {return m_quality;   }
    Double_t      getEnergy()           const {return m_energy;    }
    UInt_t        getLayer()            const {return m_firstLayer;}
    UInt_t        getTower()            const {return m_itower;    }
    TVector3      getPosition()         const {return m_position;  }
    TVector3      getDirection()        const {return m_direction; }
//    TkrParams     getTrackPar()    const;
//    double        getTrackParZ()   const {return m_position.z();}
//    TkrFitMatrix  getTrackCov()    const; 

    //! Provide info to any hit on the track
    UInt_t             getNumHits()           {return m_hits.size();   }
    TkrCandHitIter     getHitIterBegin()      {return m_hits.begin();  }
    TkrCandHitIter     getHitIterEnd()        {return m_hits.end();    }

private:
    Double_t   m_energy;
    Double_t   m_quality;

    //! Track ID 
    UInt_t     m_id;
    UInt_t     m_firstLayer;
    UInt_t     m_itower;    

    //Base information for this candidate track
    TVector3   m_position;
    TVector3   m_direction;

    //TkrCandHitVec m_hits;
    std::vector<TkrCandHit> m_hits;

    ClassDef(TkrCandTrack,1)
};
#endif
