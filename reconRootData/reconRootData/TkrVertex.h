#ifndef ROOT_TKRVERTEX_H
#define ROOT_TKRVERTEX_H

#include "TObject.h"
#include "TkrParams.h"
#include "TkrCovMat.h"
#include "TVector3.h"
#include <vector>

/** 
 * @class TkrVertex
 *
 * @brief Root object representing a reconstructed vertex
 *        Includes base summary information about the vertex
 *        and an array of track ids associated with the vertex
 *        This class should duplicate Event::TkrVertex
 * 
 * Created 13-May-2002
 * 
 * @author Tracy Usher
 *
 */

class TkrVertex : public TObject
{
public:
    TkrVertex();

    virtual ~TkrVertex();

    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    //! Define initialization method
    void initializeInfo(UInt_t id, UInt_t layer, UInt_t tower, 
        Double_t quality, Double_t energy);

    void initializeVals(const TkrParams& vtxPar, const TkrCovMat& vtxCov, 
        const TVector3& pos, const TVector3& dir);
    
    void addTrackId(UInt_t trackId) {m_trackIds.push_back(trackId);}

    //! Allow hits to be added to the list
    //void     addHit(const TkrHitPlane& hit) {m_hits.push_back(hit);}

    //! Implement methods parallel to TkrRecInfo here
    Int_t              getId()            const {return m_id;            }
    Double_t           getQuality()       const {return m_quality;       }
    Double_t           getEnergy()        const {return m_energy;        }
    Int_t              getLayer()         const {return m_firstLayer;    }
    Int_t              getTower()         const {return m_itower;        }
    const TVector3&    getPosition()      const {return m_position;      }
    const TVector3&    getDirection()     const {return m_direction;     }
    const TkrParams&   getTrackPar()      const {return m_vertexPar;     }
    Double_t           getTrackParZ()     const {return m_position.Z();  }
    const TkrCovMat&   getTrackCov()      const {return m_vertexCov;     }

    //! Provide info to any hit on the track
    UInt_t             getNumTracks() const          {return m_trackIds.size();  }
    UInt_t             getTrackId(UInt_t idx) const  {return m_trackIds[idx];    }

private:

    Double_t   m_energy;
    Double_t   m_quality;

    //! Track ID 
    Int_t     m_id;
    Int_t     m_firstLayer;
    Int_t     m_itower;        

    //! Vertex member variables
    TVector3   m_position;
    TVector3   m_direction;

    TkrParams  m_vertexPar;
    TkrCovMat  m_vertexCov;

    //! List of track IDs used by this vertex
    std::vector<UInt_t> m_trackIds;

    ClassDef(TkrVertex,1)
};
#endif
