
#ifndef TkrTrack_H
#define TkrTrack_H 

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
* $Header$
*/

#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/IInterface.h"
#include "Event/Recon/TkrRecon/TkrTrackHit.h"

// Declare Gaudi object interface ID
static const CLID& CLID_TkrTrack = InterfaceID("TkrTrack",  1, 0);

namespace Event {  // NameSpace

class TkrTrack: public TkrTrackHitVec, virtual public ContainedObject
{    
public:
    /// Constructor/destructor for the class
    TkrTrack();
   ~TkrTrack();

    //! Retrieve pointer to class defininition structure
   virtual const CLID& clID() const   { return TkrTrack::classID(); }
   static const CLID& classID()       { return CLID_TkrTrack; }

    /// Status word bits organized like:
    ///        |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         < Pat Rec Info  > <Pass > < E-Loss> < Track Energy >  <Track Fit Status>
    enum StatusBits {FOUND    = 0x0001,  //Set if track has been "found" by pat rec
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
					 BOTTOM   = 0x8000}; //Set if track traj. intercepts first Cal layer
    
    /// Utility 
    std::ostream& fillStream( std::ostream& s ) const;

    /// Access to primary quantities on track quality and scattering info
    inline unsigned int getStatusBits()          const {return m_statusBits;}
    inline double       getChiSquareFilter()     const {return m_chiSquareFilter;}
    inline double       getChiSquareSmooth()     const {return m_chiSquareSmooth;}
    inline int          getNDegreesOfFreedom()   const {return m_nDegreesOfFreedom;}
    inline double       getQuality()             const {return m_Quality;}
    ///returns the rms of the residuals between the track hits and their fitted position
    inline double       getScatter()             const {return m_rmsResid;}
    ///returns the rms of the scattering deviations of the track
    inline double       getKalThetaMS()          const {return m_KalmanThetaMS;}
    ///returns the energy estimate from Kalman Filter Fitter
    inline double       getKalEnergy()           const {return m_KalmanEnergy;}
    inline double       getKalEnergyError()      const {return m_KalmanEnergyErr;}

    /// Access to derived information on gaps
    int                 getNumGaps()             const {return m_Xgaps + m_Ygaps;}
    inline int          getNumXGaps()            const {return m_Xgaps;}
    inline int          getNumYGaps()            const {return m_Ygaps;}
    inline int          getNumXFirstGaps ()      const {return m_XistGaps;}
    inline int          getNumYFirstGaps ()      const {return m_YistGaps;}

    /// Access to fit specific information
    inline Point        getInitialPosition()     const {return m_initialPosition;}
    inline Vector       getInitialDirection()    const {return m_initialDirection;}
    inline double       getInitialEnergy()       const {return m_initialEnergy;}
    /// JCT: THE FOLLOWING SHOULD BE COMMENTED
    inline int          getNumSegmentPoints()    const {return m_numSegmentPoints;}
    inline double       chiSquareSegment(double penaltyGap = 0.)  
                                           const {return m_chisqSegment + penaltyGap*getNumGaps();}
    inline int          getNumXHits()            const {return m_nxHits;}
    inline int          getNumYHits()            const {return m_nyHits;}
	inline int          getNumFitHits()          const {return m_nxHits + m_nyHits;}
    /// JCT: THE FOLLOWING SHOULD BE COMMENTED
    inline double       getTkrCalRadlen()        const {return m_TkrCal_radlen;}

    // Access to the hit information (why must I do it this way???)
    int                 getNumHits()             const {return size();}

    /// Set functions for initializing the data members
    inline void   setInitialPosition(const Point x)   {m_initialPosition   = x;}
    inline void   setInitialDirection(const Vector x) {m_initialDirection  = x;}
    inline void   setInitialEnergy(double x)          {m_initialEnergy     = x;}
    inline void   setChiSquareFilter(double x)        {m_chiSquareFilter   = x;}
    inline void   setChiSquareSmooth(double x)        {m_chiSquareSmooth   = x;}
    inline void   setNDegreesOfFreedom(const int i)   {m_nDegreesOfFreedom = i;}
    inline void   setChiSqSegment(double x)           {m_chisqSegment      = x;}
    inline void   setQuality(double x)                {m_Quality           = x;}
    inline void   setScatter(double x)                {m_rmsResid          = x;}
    inline void   setKalThetaMS(double x)             {m_KalmanThetaMS     = x;}
    inline void   setKalEnergy(double x)              {m_KalmanEnergy      = x;}
    inline void   setKalEnergyError(double x)         {m_KalmanEnergyErr   = x;}
    inline void   setNumXGaps(int i)                  {m_Xgaps             = i;}
    inline void   setNumYGaps(int i)                  {m_Ygaps             = i;}
    inline void   setNumXFirstGaps(int i)             {m_XistGaps          = i;}
    inline void   setNumYFirstGaps(int i)             {m_YistGaps          = i;}
    inline void   setNumSegmentPoints(int i)          {m_numSegmentPoints  = i;}
    inline void   setNumXHits(int i)                  {m_nxHits            = i;}
    inline void   setNumYHits(int i)                  {m_nyHits            = i;}
    inline void   setTkrCalRadLen(double x)           {m_TkrCal_radlen     = x;}
    inline void   setStatusBit(unsigned int status)   {m_statusBits       |= status;}
	inline void   clearStatusBits()                   {m_statusBits        = 0;}
    inline void   clearEnergyStatusBits()             {m_statusBits       &= 0xff0f;}

private:	
    /// Status
    unsigned int m_statusBits;

    /// The input energy, and current position and direction
    double       m_initialEnergy;       // Initial energy at the start of the track
    Point        m_initialPosition;     // Initial position at the start of the track
    Vector       m_initialDirection;    // Initial direction at the start of the track

    // Track quality information
    double       m_chiSquareFilter;     // Chi-square/dof from Filter stage of fit
    double       m_chiSquareSmooth;     // "Smoothed" track chi-square/dof
    int          m_nDegreesOfFreedom;   // Number of degrees of freedom for above
    double       m_rmsResid;            // RMS residuals of track hits to fit
    double       m_Quality;             // Track "Quality" derived from hit counts & chisq.

    // Information from the fitter
    double       m_KalmanEnergy;        // Energy estimate from Kalman Filter Fitter
    double       m_KalmanThetaMS;       // rms scattering deviation of track

    // Hit gap information
    int          m_Xgaps;               // Number of x-meas points on track NOT used in fit
    int          m_Ygaps;               // Number of y-meas points on track Not used in fit
    int          m_XistGaps;            // Number of x-meas points in first part of track not used
    int          m_YistGaps;            // Number of y-meas points in first part of track not used

    /// Kalman Filter Track data
    int          m_numSegmentPoints;    // Effective number of 3D segments that contribute
	                                    //   to track direction
    double       m_chisqSegment;        // Chi-square for this portion of the track
    int          m_nxHits;              // Number of x meas. points USED in fit
    int          m_nyHits;              // Number of y meas. points USED in fit
    double       m_KalmanEnergyErr;     // Estimated Error on Kalman Energy
    double       m_TkrCal_radlen;       // Integrated Tracker radiation lengths 
	                                    // (uses starting track param. trajectory)
};

//typedef for the Container
typedef ObjectVector<TkrTrack>      TkrTrackCol;
typedef TkrTrackCol::iterator       TkrTrackColPtr;
typedef TkrTrackCol::const_iterator TkrTrackColConPtr;

}; //Namespace

#endif
