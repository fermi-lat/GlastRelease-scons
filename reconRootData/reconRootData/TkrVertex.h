#ifndef ROOT_TKRVERTEX_H
#define ROOT_TKRVERTEX_H

#include "TObject.h"
#include "TkrParams.h"
#include "TkrCovMat.h"
#include <vector>

#ifdef WIN32
using namespace std;
#endif

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

    virtual ~TkrVertex() {}

    //! Define initialization method
    void initializeInfo(UInt_t id, UInt_t layer, UInt_t tower, 
        Double_t quality, Double_t energy);

    void initializeVals(const TkrParams& vtxPar, const TkrCovMat& vtxCov, 
        const TVectorD& pos, const TVectorD& dir);
    
    void addTrackId(UInt_t trackId) {m_trackIds.push_back(trackId);}

    //! Allow hits to be added to the list
    //void     addHit(const TkrHitPlane& hit) {m_hits.push_back(hit);}

    //! Implement methods parallel to TkrRecInfo here
    Double_t           getQuality()       const {return m_quality;       }
    Double_t           getEnergy()        const {return m_energy;        }
    UInt_t             getLayer()         const {return m_firstLayer;    }
    UInt_t             getTower()         const {return m_itower;        }
    TVectorD           getPosition()      const {return m_position;      }
    TVectorD           getDirection()     const {return m_direction;     }
    TkrParams          getTrackPar()      const {return m_vertexPar;     }
    Double_t           getTrackParZ()     const {return m_position(2);   }
    TkrCovMat          getTrackCov()      const {return m_vertexCov;     }

    //! Provide info to any hit on the track
    UInt_t             getNumTracks()           {return m_trackIds.size();  }
    UInt_t             getTrackId(UInt_t idx)   {return m_trackIds[idx];    }

private:
    //! Track ID 
    UInt_t     m_id;

    //! Vertex member variables
    TkrParams  m_vertexPar;
    TkrCovMat  m_vertexCov;
    TVectorD   m_position;
    TVectorD   m_direction;
    Double_t   m_energy;
    Double_t   m_quality;
    UInt_t     m_firstLayer;
    UInt_t     m_itower;        

    //! List of track IDs used by this vertex
    std::vector<UInt_t> m_trackIds;

    ClassDef(TkrVertex,1)
};
#endif
