
#ifndef __TkrKalFitTrack_H
#define __TkrKalFitTrack_H 1

#include <vector>
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrFitTrackBase.h"

/** 
* @class TkrKalFitTrack
*
* @brief Contains all necessary information to define a fit track in GLAST
*        The class contains basic information resulting from the fit of the 
*        track (track quality, RMS residuals, scattering, gaps, etc.) as 
*        well as a vector of "hits." Each "hit" contains fit information for
*        each plane with a measurement. 
*
* Adapted from original version of Bill Atwood
*
* @author The Tracking Software Group
*
* $Header$
*/
namespace Event {  // NameSpace

class TkrKalFitTrack: public TkrFitTrackBase
{    
public:
    /// Constructor/destructor for the class
    TkrKalFitTrack();
   ~TkrKalFitTrack();

    enum          Status {EMPTY, FOUND, CRACK}; 

    void     initializeBase(Status stat, int type, double energy0, Point& x0, Vector& dir);
    void     initializeQual(double chiSq, double ChiSqSmooth, double rms, double quality, 
                            double e, double e_err, double ms);
    void     initializeGaps(int xgaps, int ygaps, int x1st, int y1st);
    void     initializeKal( int nSegPoints, double chisqSeg, int nxHits, int nyHits, 
                            double radlen);

    /// Define the TkrRecInfo access methods
    /// Provides access to the basic information needed by external users
    double        getQuality()                       const {return m_Q;};
    double        getEnergy(TrackEnd end = Start)    const {
        return getFoLPlane(end).getEnergy();
    }
    int           getLayer(TrackEnd end = Start)     const {
        return getFoLPlane(end).getIDPlane();
    }
    int           getTower(TrackEnd end = Start)     const {
        return getFoLPlane(end).getIDTower();
    }
    Point         getPosition(TrackEnd end = Start)  const {
        return getFoLPlane(end).getPoint(TkrFitHit::SMOOTH);
    }

    Vector        getDirection(TrackEnd end = Start) const;
    Ray           getRay(TrackEnd end = Start)       const;

    TkrFitPar     getTrackPar(TrackEnd end = Start)  const {
        return getFoLPlane(end).getHit(TkrFitHit::SMOOTH).getPar();
    }
    double        getTrackParZ(TrackEnd end = Start) const {
        return getFoLPlane(end).getZPlane();
    }
    TkrFitMatrix  getTrackCov(TrackEnd end = Start)  const {
        return getFoLPlane(end).getHit(TkrFitHit::SMOOTH).getCov();
    }
    bool          empty(int numHits)                 const;
    
    /// Utility 
    void   writeOut(MsgStream& log) const; 

    /// Access to primary quantities on track quality and scattering info
    inline Status status()                 const {return m_status;}
    inline double getChiSquare()           const {return m_chisq;}
    inline double getChiSquareSmooth()     const {return m_chisqSmooth;}
    inline double getScatter()             const {return m_rmsResid;}
    inline double getKalThetaMS()          const {return m_KalThetaMS;}
    inline double getKalEnergy()           const {return m_KalEnergy;}
    inline double getKalEnergyError()      const {return m_KalEnergyErr;}

    /// Access to derived information on gaps
    int           getNumGaps()             const {return m_Xgaps + m_Ygaps;}
    inline int    getNumXGaps()            const {return m_Xgaps;}
    inline int    getNumYGaps()            const {return m_Ygaps;}
    inline int    getNumXFirstGaps ()      const {return m_XistGaps;}
    inline int    getNumYFirstGaps ()      const {return m_YistGaps;}

    /// Access to fit specific information
    inline Point  getInitialPosition()     const {return m_x0;}
    inline Vector getInitialDirection()    const {return m_dir;}
    inline Point  getPosAtZ(double deltaZ) const {return m_x0+deltaZ*m_dir;} 
    inline double getStartEnergy()         const {return m_energy0;}
    inline int    getNumSegmentPoints()    const {return m_numSegmentPoints;}
    inline double chiSquareSegment(double penaltyGap = 0.)  
                                           const {return m_chisqSegment + penaltyGap*getNumGaps();}
    inline int    getNumXHits()            const {return m_nxHits;}
    inline int    getNumYHits()            const {return m_nyHits;}
    inline int    getType()                const {return m_type;}
    inline double getTkrCalRadlen()        const {return m_TkrCal_radlen;}

    /// Access to derived information on kinks
    double        getKink(int iplane)      const;
    double        getKinkNorma(int iplane) const;

    /// Access to plane information
    int           getNumHits()             const {return size();}
    TkrFitPlane   getFoLPlane(TrackEnd end = Start) const;

    /// Set functions for initializing the data members
    inline void   setInitialPosition(const Point x)   {m_x0 = x;}
    inline void   setInitialDirection(const Vector x) {m_dir = x;}
    inline void   setChiSquare(double x)              {m_chisq = x;}
    inline void   setChiSquareSmooth(double x)        {m_chisqSmooth = x;}
    inline void   setChiSqSegment(double x)           {m_chisqSegment = x;}
    inline void   setQuality(double x)                {m_Q = x;}
    inline void   setScatter(double x)                {m_rmsResid = x;}
    inline void   setKalThetaMS(double x)             {m_KalThetaMS = x;}
    inline void   setKalEnergy(double x)              {m_KalEnergy = x;}
    inline void   setKalEnergyError(double x)         {m_KalEnergyErr = x;}
    inline void   setNumXGaps(int i)                  {m_Xgaps = i;}
    inline void   setNumYGaps(int i)                  {m_Ygaps = i;}
    inline void   setNumXFirstGaps(int i)             {m_XistGaps = i;}
    inline void   setNumYFirstGaps(int i)             {m_YistGaps = i;}
    inline void   setStartEnergy(double x)            {m_energy0 = x;}
    inline void   setNumSegmentPoints(int i)          {m_numSegmentPoints = i;}
    inline void   setNumXHits(int i)                  {m_nxHits = i;}
    inline void   setNumYHits(int i)                  {m_nyHits = i;}
    inline void   setType(int i)                      {m_type = i;}
    inline void   setTkrCalRadLen(double x)           {m_TkrCal_radlen = x;}
    inline void   setStatus(Status status)            {m_status = status;}
private:	
    /// Status
    Status m_status;
    int    m_type; 

    /// The input energy, and current position and direction
    double m_energy0;
    Point  m_x0;
    Vector m_dir;

    // Track quality information
    double m_chisq;                // Standard track chi-square
    double m_chisqSmooth;          // "Smoothed" track chi-square
    double m_rmsResid;             // RMS residuals of track hits to fit
    double m_Q;                    // Track "Quality"

    // Information from the fitter
    double m_KalEnergy;
    double m_KalThetaMS;

    // Hit gap information
    int    m_Xgaps;
    int    m_Ygaps;
    int    m_XistGaps;
    int    m_YistGaps;

    /// KalTrack data
    int    m_numSegmentPoints;
    double m_chisqSegment;
    int    m_nxHits;
    int    m_nyHits;
    double m_KalEnergyErr;
    double m_TkrCal_radlen; 
};

}; //Namespace

#endif
