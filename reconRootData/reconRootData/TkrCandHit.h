#ifndef ROOT_TKRCANDHIT_H
#define ROOT_TKRCANDHIT_H

#include "TObject.h"
#include "TVector3.h"
#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif

/** 
 * @class TkrCandHit
 * @brief Root object containing track candidate hit information
 *
 * This class should exactly duplicate the Event::TkrPatCandHit TDS class.
 * 
 * @author Tracy Usher
 *
 * $Header$
 */

class TkrCandHit : public TObject
{
public:
    //! This probably does not want to be 
    enum AXIS{X,Y,NONE};

    TkrCandHit();

    TkrCandHit(const TVector3& pos, UInt_t hitId, UInt_t tower, UInt_t layer, AXIS view);

    virtual ~TkrCandHit() {}

    void Clear(Option_t *option="");
     
    void Print(Option_t *option="") const;

    //! Provide an initialize method (in the case of the default constructor)
    void     initialize(const TVector3& pos, UInt_t hitId, UInt_t tower, UInt_t layer, AXIS view);

    //! Group the methods for returning info here
    const TVector3&  getPosition()  const {return m_position; }
    UInt_t           getHitIndex()  const {return m_hitIndex; }
    UInt_t           getTower()     const {return m_tower;    }
    UInt_t           getPlane()     const {return m_plane;    }
    AXIS             getView()      const {return m_view;     }

private:
    UInt_t   m_hitIndex;          //Index into TkrSiClusters list of this hit
    UInt_t   m_tower;             //Index of the tower containing this hit
    UInt_t   m_plane;             //Plane (layer) - TkrRecon convention
    AXIS     m_view;              //Primary measurement projection (x or y)};
    TVector3 m_position;          //xyz coordinates of this hit

    ClassDef(TkrCandHit,1)
};

typedef std::vector<TkrCandHit>                 TkrCandHitVec;
typedef std::vector<TkrCandHit>::const_iterator TkrCandHitIter;
#endif
