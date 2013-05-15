#ifndef ROOT_TKRTRACKHIT_H
#define ROOT_TKRTRACKHIT_H

#include "TObject.h"
#include "commonRootData/idents/TkrId.h"
#include "TkrCluster.h"
#include "TkrTrackParams.h"
#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif

/** 
 * @class TkrTrackHit
 * @brief Root object containing fit track information for one plane
 * 
 * It contains summary information for the layer plus track parameters
 * and covariance matrix for each type of hit - measured, predicted, fit (filtered),
 * and smoothed. It also keeps the scattering matrix between planes.
 * This class should duplicate Event::TkrFitPlane TDS class.
 * 
 * @author Tracy Usher
 *
 */

class TkrTrackHit : public TObject
{
public:

  /** @enum The ParamType enum defines the possible "types" of track parameters:
   * - MEASURED  = Measured values from the tracker Si planes
   * - PREDICTED = Projected estimated values from previous plane
   * - FILTERED  = Parameters from track fit (KF filter stage)
   * - SMOOTHED  = Parameters from Kalman Filter smoothing stage
   * - QMATERIAL = For access to the contribution from scattering
   * - UNKNOWN}  = Unknown
   */
    enum ParamType  {MEASURED,               // Measured values from the tracker Si planes
                     PREDICTED,              // Projected estimated values from previous plane
                     FILTERED,               // Parameters from track fit (KF filter stage)
                     REVFIT,                 // Parameters from track fit (KF filter stage)
                     SMOOTHED,               // Parameters from Kalman Filter smoothing stage
                     QMATERIAL,              // For access to the contribution from scattering
                     UNKNOWN};               // Unknown

    /// Status word bits organized like:
    /// low:   |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         < volume info  >  <   hit type    > <    track fitting status          >
    /// high:  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///                                                                < more hit type >
    enum StatusBits {HITONFIT     = 0x0001,  // Hit is used in the fit
                     HASMEASURED  = 0x0002,  // Hit has valid measured parameters
                     HASPREDICTED = 0x0004,  // Hit has valid predicted parameters
                     HASFILTERED  = 0x0008,  // Hit has valid filtered parameters
                     HASSMOOTHED  = 0x0010,  // Hit has valid smoothed parameters
                     HASMATERIAL  = 0x0020,  // Hit has valid material matrix
                     UPWARDS      = 0x0040,  // Track direction is upwards (tz > 0)

                     HITISSSD     = 0x0100,  // Hit comes from a SSD
                     HITISDEADST  = 0x0200,  // Hit coresponds to a dead SSD Strip
                     HITISGAP     = 0x0400,  // Hit comes from gap between SSDs
                     HITISTWR     = 0x0800,  // Hit comes outside live SSD plane

                     MEASURESX    = 0x1000,  // Plane measures in X direction
                     MEASURESY    = 0x2000,  // Plane measures in Y direction
                     HASVALIDTKR  = 0x8000,  // Valid track volume identifier

                     HITISUNKNOWN = 0x10000, // Missing cluster, but fails all tests
                     HITISDEADPLN = 0x20000  // Entire plane is dead
    };

    /// Default (null) constructor (just in case...)
    TkrTrackHit();

    /// Construct all but the track parameters, they must be set during recon stage
    TkrTrackHit(TkrCluster*           cluster,
                commonRootData::TkrId tkrID,
                const Double_t          hitZ,
                const Double_t          hitEnergy,
                const Double_t          hitRadLen,
                const Double_t          hitActDist,
                const Double_t          hitChiFilter,
                const Double_t          hitChiSmooth) :
                    m_statusBits(0), m_cluster(cluster), m_hitID(tkrID), m_zPlane(hitZ), 
                    m_energy(hitEnergy), m_radLen(hitRadLen), m_activeDist(hitActDist),
                    m_chiSquareFilter(hitChiFilter), m_chiSquareSmooth(hitChiSmooth)
    {
        if (m_cluster != 0)  m_statusBits  = HITISSSD | HITONFIT;
        if (tkrID.hasTray()) m_statusBits |= HASVALIDTKR;
    }

    virtual ~TkrTrackHit() {}

    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TkrTrackHit &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    //! Two initialize methods needed for this class
    void     initializeInfo(TkrCluster*           cluster,
                            commonRootData::TkrId tkrID,
                            const Double_t        hitZ,
                            const Double_t        hitEnergy,
                            const Double_t        hitRadLen,
                            const Double_t        hitActDist,
                            const Double_t        hitChiFilter,
                            const Double_t        hitChiSmooth);

    void     initializeHits(const TkrTrackParams& meas, const TkrTrackParams& pred, 
                            const TkrTrackParams& fit,  const TkrTrackParams& smooth, 
                            const TkrTrackParams& material);

    /// Access the status bits to determine details of the hit
    inline const UInt_t getStatusBits() const {return m_statusBits;}

    /// Answer quick questions based on status bits
    inline const bool   validCluster()      const {return  m_cluster != 0;}
        inline const bool   hitUsedOnFit()      const {return (m_statusBits & HITONFIT    ) == HITONFIT;}
    inline const bool   validMeasuredHit()  const {return (m_statusBits & HASMEASURED ) == HASMEASURED;}
    inline const bool   validPredictedHit() const {return (m_statusBits & HASPREDICTED) == HASPREDICTED;}
    inline const bool   validFilteredHit()  const {return (m_statusBits & HASFILTERED ) == HASFILTERED;}
    inline const bool   validSmoothedHit()  const {return (m_statusBits & HASSMOOTHED ) == HASSMOOTHED;}
    inline const bool   validMaterial()     const {return (m_statusBits & HASMATERIAL ) == HASMATERIAL;}
    inline const bool   upwards()           const {return (m_statusBits & UPWARDS )     == UPWARDS;}

    /// Access the hit's associated information
    inline const TkrCluster* getClusterPtr()      const {return m_cluster;   }
    inline const Double_t    getZPlane()          const {return m_zPlane;    }
    inline const Double_t    getEnergy()          const {return m_energy;    }
    inline const Double_t    getRadLen()          const {return m_radLen;    }
    inline const Double_t    getActiveDist()      const {return m_activeDist;}
    inline const Double_t    getChiSquareFilter() const {return m_chiSquareFilter;}
    inline const Double_t    getChiSquareRevFit() const {return m_chiSquareRevFit;}
    inline const Double_t    getChiSquareSmooth() const {return m_chiSquareSmooth;}
    inline const commonRootData::TkrId getTkrId() const {return m_hitID;          }

    /// @enum Allow rudimentary access to the hit information here 
    /// Method to return index to give TkrTrackParams for returning specifc parameter
    enum SSDDirection  {SSDMEASURED, SSDNONMEASURED};

    const int getParamIndex(TkrTrackHit::SSDDirection meas, TkrTrackParams::ParamType type) const; 
    int       getParamIndex(TkrTrackHit::SSDDirection meas, TkrTrackParams::ParamType type); 

    const TVector3           getPoint(TkrTrackHit::ParamType type)      const;
    TVector3                 getPoint(TkrTrackHit::ParamType type);
    const TVector3           getDirection(TkrTrackHit::ParamType type)  const;
    TVector3                 getDirection(TkrTrackHit::ParamType type);
    const Double_t           getMeasuredPosition(TkrTrackHit::ParamType type) const;
    const Double_t           getMeasuredSlope(TkrTrackHit::ParamType type) const;
    const Double_t           getNonMeasuredPosition(TkrTrackHit::ParamType type) const;
    const Double_t           getNonMeasuredSlope(TkrTrackHit::ParamType type) const;

    /// Direct access to track params
    const TkrTrackParams& getTrackParams(ParamType type) const;
    TkrTrackParams& getTrackParams(ParamType type);

    /// Set methods used during pattern recognition/track fit stages
    inline void setClusterPtr(TkrCluster* cluster)     {m_cluster         = cluster;}
    inline void setTkrID(commonRootData::TkrId& tkrID) {m_hitID           = tkrID;  }
    inline void setEnergy(const Double_t energy)       {m_energy          = energy; }
    inline void setZPlane(const Double_t z)            {m_zPlane          = z;}
    inline void setRadLen(const Double_t rl)           {m_radLen          = rl;}
    inline void setActiveDist(const Double_t d)        {m_activeDist      = d;}
    inline void setChiSquareFilter(const Double_t c)   {m_chiSquareFilter = c;}
    inline void setChiSquareRevFit(const Double_t c)   {m_chiSquareRevFit = c;}
    inline void setChiSquareSmooth(const Double_t c)   {m_chiSquareSmooth = c;}

    inline void setStatusBit(UInt_t bitToSet)          {m_statusBits |=  bitToSet;}
    inline void clearStatusBit(UInt_t bitToClear)      {m_statusBits &= ~bitToClear;}

    /// Override the "new" operator in order to use the ReconObjectManager pool management
    void* operator new(size_t size);

    void* operator new(size_t size, void* vp);

    /// Override the "delete" operator in order to use the ReconObjectManager pool management
    void  operator delete(void* p);

private:
    inline const Double_t getCoordinate(const TkrTrackParams& params, int coord) const;
    inline Double_t       getCoordinate(const TkrTrackParams& params, int coord);

    UInt_t                m_statusBits;      // See StatusBits enumeration above for definitions
    TkrCluster*           m_cluster;         // Pointer to the cluster associated with this hit
    commonRootData::TkrId m_hitID;           // Complete TkrId identifying the details of this hit/plane
    Double32_t            m_zPlane;          // Z location of plane
    Double32_t            m_energy;          // Energy of track at this plane
    Double32_t            m_radLen;          // Radiation Lengths encountered from the previous hit
    Double32_t            m_activeDist;      // The distance inside (positive) hit SSD (neg. if outside)
    Double32_t            m_chiSquareFilter; // hit chi-square at filter stage of fit
    Double32_t            m_chiSquareRevFit; // hit chi-square at filter stage of fit
    Double32_t            m_chiSquareSmooth; // hit chi-square at smooth stage of fit
        
    TkrTrackParams        m_hitMeas;
    TkrTrackParams        m_hitPred;
    TkrTrackParams        m_hitFit;
    TkrTrackParams        m_hitRevFit;
    TkrTrackParams        m_hitSmooth;
    TkrTrackParams        m_Qmaterial;  // holds the covariance matrix of the material effects 

    ClassDef(TkrTrackHit,3)
};

//typedef std::vector<TkrTrackHit>                 TkrTrackHitVector;
//typedef std::vector<TkrTrackHit>::const_iterator TkrTrackHitIter;
#endif
