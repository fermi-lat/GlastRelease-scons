#ifndef ROOT_TKRHITPLANE_H
#define ROOT_TKRHITPLANE_H

#include "TObject.h"
#include "TkrFitHit.h"
#include "TkrCovMat.h"
#include <vector>
#ifdef WIN32
using namespace std;
#endif

/** 
 * @class TkrHitPlane
 *
 * @brief Root object containing fit track information for one plane
 *        It contains summary information for the layer plus track parameters
 *        and covariance matrix for each type of hit - measured, predicted, fit (filtered),
 *        and smoothed. It also keeps the scattering matrix between planes.
 *        This class should duplicate Event::TkrFitPlane
 * 
 * Created 13-May-2002
 * 
 * @author Tracy Usher
 *
 */

class TkrHitPlane : public TObject
{
public:
    //! This probably does not want to be 
    enum AXIS{X,Y,NONE};

    TkrHitPlane();

    virtual ~TkrHitPlane() {}

    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    //! Two initialize methods needed for this class
    void     initializeInfo(UInt_t hit, UInt_t tower, UInt_t plane, AXIS proj,
                            AXIS nextProj, Double_t z, Double_t energy, Double_t radLen, Double_t activeDist);
    void     initializeHits(const TkrFitHit& meas, const TkrFitHit& pred, const TkrFitHit& fit,
                            const TkrFitHit& smooth, const TkrCovMat& material);

    //! Methods for retrieving the information
    UInt_t           getIdHit()              const {return m_IdHit;         }
    UInt_t           getIdTower()            const {return m_IdTower;       }
    UInt_t           getIdPlane()            const {return m_IdPlane;       }
    AXIS             getProjection()         const {return m_Projection;    }
    AXIS             getProjPlus()           const {return m_ProjPlus;      }
    Double_t         getZplane()             const {return m_Zplane;        }
    Double_t         getEnePlane()           const {return m_EnePlane;      }
    Double_t         getRadLen()             const {return m_RadLen;        }
    Double_t         getActiveDist()         const {return m_activeDist;    }

    const TkrFitHit& getHitMeas()            const {return m_HitMeas;       }
    const TkrFitHit& getHitPred()            const {return m_HitPred;       }
    const TkrFitHit& getHitFit()             const {return m_HitFit;        }
    const TkrFitHit& getHitSmooth()          const {return m_HitSmooth;     }
    const TkrCovMat& getQmaterial()          const {return m_Qmaterial;     }

private:
    Double_t   m_Zplane;      // Z coordinate of this plane
    Double_t   m_EnePlane;    // Fitter energy at this plane
    Double_t   m_RadLen;      // Radiation lengths to this plane
    Double_t   m_activeDist;  // Distance to edge of nearest active volume

    UInt_t     m_IdHit;       // Index of cluster hit in the TkrSiClusters 
    UInt_t     m_IdTower;     // Tower ID
    UInt_t     m_IdPlane;     // Plane ID (TkrRecon convention)
    AXIS       m_Projection;  // X or Y measuring SSD plane
    AXIS       m_ProjPlus;    // same but for "next" plane

    TkrFitHit  m_HitMeas;     // Measured track parameters
    TkrFitHit  m_HitPred;     // Predicted track parameters
    TkrFitHit  m_HitFit;      // "Fit" track parameters
    TkrFitHit  m_HitSmooth;   // "Smoothed" track parameters

    TkrCovMat  m_Qmaterial;   // Covariance matrix due to multiple scatter

    ClassDef(TkrHitPlane,2)
};

typedef std::vector<TkrHitPlane>                 TkrHitPlaneVector;
typedef std::vector<TkrHitPlane>::const_iterator TkrHitPlaneIter;
#endif
