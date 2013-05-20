#ifndef ROOT_TKRTRACK_H
#define ROOT_TKRTRACK_H

#include "TObjArray.h"
#include "TkrTrackHit.h"
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif

/** 
* @class TkrTrack
*
* @brief TDS Data Object defining a Tracker reconstructred track
*
* This class defines the basic TkrRecon structure containing output information
* for a track reconstructed in the LAT Tracker. 
* This information includes:
* 1) Reconstruction Status information (see StatusBits)
* 2) Summary reconstruction information including:
*    a) Filter and Smooth step chi-squares
*    b) Track ordering information (primarily from pat rec)
*    c) Track rms, estimated energy, ...
*    d) Variables to help estimate track quality
*    e) etc.
* 3) An ObjectVector containing information for all hits on the track
*    (TkrTrackHit), including the resulting track parameters and 
*    covariance matrices (TkrTrackParams)
*
* @author The Tracking Software Group
*
 * 
 * $Header$
 */

class TkrTrack : public TObjArray
{
public:

    /// @enum Status word bits organized like:
    ///        |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         < Pat Rec Info  > <Pass > < E-Loss> < Track Energy >  <Track Fit Status>
    enum StatusBits {
                    FOUND    = 0x0001,  //Set if track has been "found" by pat rec
                    FILTERED = 0x0002,  //Set if track fit filter stage has been run
                    SMOOTHED = 0x0004,  //Set if track fit smoother has been run
                    REVFILTR = 0x0008,  //Set if track has been reverse-filtered
                    CALENERGY= 0x0010,  //Set if track energy from raw calorimeter info
                    LATENERGY= 0x0020,  //Set if track energy from TKR+CAL constrained
                    USERENERGY= 0x0040, //Set if track energy set by user
                    MCENERGY = 0x0080,  //Set if energy from users or from MC truth
                    RADELOSS = 0x0100,  //Set if radiative energy loss used (e+/e- fitting)
                    MIPELOSS = 0x0200,  //Set if Bethe-Block energy loss used (not e+/e-)
                    ONEPASS  = 0x0400,  //Set if the full first pass track fit finished
                    TWOPASS  = 0x0800,  //Set if an iteration of the first fit finished
                    PRCALSRCH= 0x1000,  //Set if Pat. Rec. used Cal Energy Centroid
                    PRBLNSRCH= 0x2000,  //Set if Pat. Rec. used only Track info.
                    TOP      = 0x4000,  //Set if track traj. intercepts top tracker plane
                    BOTTOM   = 0x8000,  //Set if track traj. intercepts first Cal layer

                    GHOST      = 0x010000, // set if track contains ghost clusters
                    DIAGNOSTIC = 0x020000, // set if track contains diagnostic ghost clusters
                    COSMICRAY  = 0x040000, // set if the track was found by the cosmic-ray finder (don't use these for gamma-ray analysis!)
                    TREEBASED  = 0x080000, // set if the track was found using the Tree Based pattern recognition

                    COMPOSITE  = 0x100000, // set if "hits" are composite positions instead of clusters
                    HASKINKS   = 0x200000  // set if the track has kinks in it
    };

    TkrTrack();

    virtual ~TkrTrack();

    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TkrTrack &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    //! Define two initialization methods
    //void     initializeInfo(UInt_t id, UInt_t xgaps, UInt_t ygaps, 
    //    UInt_t x1st, UInt_t y1st);
    //void     initializeQual(Double32_t chiSq, Double32_t ChiSqSmooth, 
    //    Double32_t rms, Double32_t quality, Double32_t e, Double32_t ms);

    /// Access to primary quantities on track quality and scattering info
    inline UInt_t       getStatusBits()          const {return m_statusBits;}
    inline Double_t     getChiSquareFilter()     const {return m_chiSquareFilter;}
    inline Double_t     getChiSquareSmooth()     const {return m_chiSquareSmooth;}
    inline UInt_t       getNDegreesOfFreedom()   const {return m_nDegreesOfFreedom;}
    inline Double_t     getQuality()             const {return m_Quality;}
    ///returns the rms of the residuals between the track hits and their fitted position
    inline Double_t     getScatter()             const {return m_rmsResid;}
    ///returns the rms of the scattering deviations of the track
    inline Double_t     getKalThetaMS()          const {return m_KalmanThetaMS;}
    ///returns the energy estimate from Kalman Filter Fitter
    inline Double_t     getKalEnergy()           const {return m_KalmanEnergy;}
    inline Double_t     getKalEnergyError()      const {return m_KalmanEnergyErr;}

    /// Access to derived information on gaps
    UInt_t              getNumGaps()             const {return m_Xgaps + m_Ygaps;}
    inline UInt_t       getNumXGaps()            const {return m_Xgaps;}
    inline UInt_t       getNumYGaps()            const {return m_Ygaps;}
    inline UInt_t       getNumXFirstGaps ()      const {return m_XistGaps;}
    inline UInt_t       getNumYFirstGaps ()      const {return m_YistGaps;}

    /// Access to fit specific information
    inline TVector3     getInitialPosition()     const {return m_initialPosition;}
    inline TVector3     getInitialDirection()    const {return m_initialDirection;}
    inline Double_t     getInitialEnergy()       const {return m_initialEnergy;}
    /// JCT: THE FOLLOWING SHOULD BE COMMENTED
    inline UInt_t       getNumSegmentPoints()    const {return m_numSegmentPoints;}
    inline Double_t     chiSquareSegment(Double_t penaltyGap = 0.)  
                                           const {return m_chisqSegment + penaltyGap*getNumGaps();}
    inline UInt_t       getNumXHits()            const {return m_nxHits;}
    inline UInt_t       getNumYHits()            const {return m_nyHits;}
        inline UInt_t       getNumFitHits()          const {return m_nxHits + m_nyHits;}
    /// JCT: THE FOLLOWING SHOULD BE COMMENTED
    inline Double_t     getTkrCalRadlen()        const {return m_TkrCal_radlen;}

    /// Access to the hits
    TIterator*          Iterator()               const {return MakeIterator();}
    UInt_t              Size()                   const {return GetEntries();}

    /// Set functions for initializing the data members
    inline void   setInitialPosition(const TVector3& x)  {m_initialPosition   = x;}
    inline void   setInitialDirection(const TVector3& x) {m_initialDirection  = x;}
    inline void   setInitialEnergy(Double_t x)           {m_initialEnergy     = x;}
    inline void   setChiSquareFilter(Double_t x)         {m_chiSquareFilter   = x;}
    inline void   setChiSquareSmooth(Double_t x)         {m_chiSquareSmooth   = x;}
    inline void   setNDegreesOfFreedom(const UInt_t i)   {m_nDegreesOfFreedom = i;}
    inline void   setChiSqSegment(Double_t x)            {m_chisqSegment      = x;}
    inline void   setQuality(Double_t x)                 {m_Quality           = x;}
    inline void   setScatter(Double_t x)                 {m_rmsResid          = x;}
    inline void   setKalThetaMS(Double_t x)              {m_KalmanThetaMS     = x;}
    inline void   setKalEnergy(Double_t x)               {m_KalmanEnergy      = x;}
    inline void   setKalEnergyError(Double_t x)          {m_KalmanEnergyErr   = x;}
    inline void   setNumXGaps(UInt_t i)                  {m_Xgaps             = i;}
    inline void   setNumYGaps(UInt_t i)                  {m_Ygaps             = i;}
    inline void   setNumXFirstGaps(UInt_t i)             {m_XistGaps          = i;}
    inline void   setNumYFirstGaps(UInt_t i)             {m_YistGaps          = i;}
    inline void   setNumSegmentPoints(UInt_t i)          {m_numSegmentPoints  = i;}
    inline void   setNumXHits(UInt_t i)                  {m_nxHits            = i;}
    inline void   setNumYHits(UInt_t i)                  {m_nyHits            = i;}
    inline void   setTkrCalRadLen(Double_t x)            {m_TkrCal_radlen     = x;}
    inline void   setStatusBit(UInt_t status)            {m_statusBits       |= status;}
        inline void   clearStatusBits()                      {m_statusBits        = 0;}
    inline void   clearEnergyStatusBits()                {m_statusBits       &= 0xff0f;}

    TkrTrackHit*        getHit(UInt_t iHit)      const {return iHit < Size() ? dynamic_cast<TkrTrackHit*>(At(iHit)) : 0; }

    /// Override the "new" operator in order to use the ReconObjectManager pool management
    void* operator new(size_t size);

    void* operator new(size_t size, void* vp);

    /// Override the "delete" operator in order to use the ReconObjectManager pool management
    void  operator delete(void* p);
    
    inline void        setRangeEnergy(Double_t x)  {m_rangeEnergy = x;}
    inline Double_t  getRangeEnergy()              const {return m_rangeEnergy;}

private:
    /// Status
    UInt_t       m_statusBits;

    /// The input energy, and current position and direction
    Double32_t   m_initialEnergy;       // Initial energy at the start of the track
    TVector3     m_initialPosition;     // Initial position at the start of the track
    TVector3     m_initialDirection;    // Initial direction at the start of the track

    // Track quality information
    Double32_t   m_chiSquareFilter;     // Chi-square/dof from Filter stage of fit
    Double32_t   m_chiSquareSmooth;     // "Smoothed" track chi-square/dof
    UInt_t       m_nDegreesOfFreedom;   // Number of degrees of freedom for above
    Double32_t   m_rmsResid;            // RMS residuals of track hits to fit
    Double32_t   m_Quality;             // Track "Quality" derived from hit counts & chisq.

    // Information from the fitter
    Double32_t   m_KalmanEnergy;        // Energy estimate from Kalman Filter Fitter
    Double32_t   m_KalmanThetaMS;       // rms scattering deviation of track

    // Hit gap information
    UInt_t       m_Xgaps;               // Number of x-meas points on track NOT used in fit
    UInt_t       m_Ygaps;               // Number of y-meas points on track Not used in fit
    UInt_t       m_XistGaps;            // Number of x-meas points in first part of track not used
    UInt_t       m_YistGaps;            // Number of y-meas points in first part of track not used

    /// Kalman Filter Track data
    UInt_t       m_numSegmentPoints;    // Effective number of 3D segments that contribute
                                            //   to track direction
    Double32_t   m_chisqSegment;        // Chi-square for this portion of the track
    UInt_t       m_nxHits;              // Number of x meas. points USED in fit
    UInt_t       m_nyHits;              // Number of y meas. points USED in fit
    Double32_t   m_KalmanEnergyErr;     // Estimated Error on Kalman Energy
    Double32_t   m_TkrCal_radlen;       // Integrated Tracker radiation lengths 
    Double32_t   m_rangeEnergy;           // energy estimated from the range 

    ClassDef(TkrTrack,6)

};
#endif
