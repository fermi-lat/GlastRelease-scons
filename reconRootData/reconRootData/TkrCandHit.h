#ifndef ROOT_TKRCANDHIT_H
#define ROOT_TKRCANDHIT_H

#include "TObject.h"
#include "TVectorD.h"
#include <vector>

/** 
 * @class TkrCandHit
 *
 * @brief Root object containing track candidate hit information
 *        This class should exactly duplicate the TkrPatCandHit class
 * 
 * Created 15-May-2002
 * 
 * @author Tracy Usher
 *
 */

class TkrCandHit : public TObject
{
public:
    //! This probably does not want to be 
    enum AXIS{X,Y,NONE};

    TkrCandHit();

    TkrCandhit(const TVectorD& pos, UInt_t hitId, UInt_t tower, UInt_t layer, AXIS view);

    virtual ~TkrCandHit() {}

    //! Provide an initialize method (in the case of the default constructor)
    void     initialize(const TVectorD& pos, UInt_t hitId, UInt_t tower, UInt_t layer, AXIS view);

    //! Group the methods for returning info here
    TVectorD         getPosition()  const {return m_position; }
    UInt_t           getHitIndex()  const {return m_hitIndex; }
    UInt_t           getTower()     const {return m_tower;    }
    UInt_t           getPlane()     const {return m_plane;    }
    AXIS             getView()      const {return m_view;     }

private:
    TVectorD m_position;          //xyz coordinates of this hit
    UInt_t   m_hitIndex;          //Index into TkrSiClusters list of this hit
    UInt_t   m_tower;             //Index of the tower containing this hit
    UInt_t   m_plane;             //Plane (layer) - TkrRecon convention
    AXIS     m_view;              //Primary measurement projection (x or y)};

    ClassDef(TkrCandHit,1)
};

typedef std::vector<TkrCandHit>                 TkrCandHitVec;
typedef std::vector<TkrCandHit>::const_iterator TkrCandHitIter;
#endif
