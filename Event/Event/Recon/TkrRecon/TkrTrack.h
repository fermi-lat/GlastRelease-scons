
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

    enum StatusBits {Found    = 0x0001,  // Set if track has been "found" by pat rec
                     Filtered = 0x0002,  // Set if track fit filter stage has been run
                     Smoothed = 0x0004,  // Set if track fit smoother has been run
                     OnePass  = 0x0100,  // Set if the full first pass track fit finished
                     TwoPass  = 0x0200}; // Set if an iteration of the first fit finished
    
    /// Utility 
    void   writeOut(MsgStream& log) const; 

    /// Access to primary quantities on track quality and scattering info
    inline unsigned int getStatusBits()          const {return m_statusBits;}
    inline double       getChiSquareFilter()     const {return m_chiSquareFilter;}
    inline double       getChiSquareSmooth()     const {return m_chiSquareSmooth;}
    inline int          getNDegreesOfFreedom()   const {return m_nDegreesOfFreedom;}
    inline double       getQuality()             const {return m_Quality;}
    inline double       getScatter()             const {return m_rmsResid;}
    inline double       getKalThetaMS()          const {return m_KalmanThetaMS;}
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
    inline int          getNumSegmentPoints()    const {return m_numSegmentPoints;}
    inline double       chiSquareSegment(double penaltyGap = 0.)  
                                           const {return m_chisqSegment + penaltyGap*getNumGaps();}
    inline int          getNumXHits()            const {return m_nxHits;}
    inline int          getNumYHits()            const {return m_nyHits;}
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
    double       m_Quality;             // Track "Quality"

    // Information from the fitter
    double       m_KalmanEnergy;        // Energy estimate from Kalman Filter Fitter
    double       m_KalmanThetaMS;       // rms scattering deviation of track

    // Hit gap information
    int          m_Xgaps;               //
    int          m_Ygaps;               //
    int          m_XistGaps;            //
    int          m_YistGaps;            //

    /// Kalman Filter Track data
    int          m_numSegmentPoints;    //
    double       m_chisqSegment;        //
    int          m_nxHits;              //
    int          m_nyHits;              //
    double       m_KalmanEnergyErr;     //
    double       m_TkrCal_radlen;       //
};

//typedef for the Container
typedef ObjectVector<TkrTrack>      TkrTrackCol;
typedef TkrTrackCol::iterator       TkrTrackColPtr;
typedef TkrTrackCol::const_iterator TkrTrackColConPtr;

}; //Namespace

#endif
