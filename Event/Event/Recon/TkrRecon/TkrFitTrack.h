
#ifndef __TkrFitTrack_H
#define __TkrFitTrack_H 1

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "Event/Recon/TkrRecon/TkrRecInfo.h"
#include "Event/Recon/TkrRecon/TkrFitPlane.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"

/** 
* @class TkrFitTrack
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

  class TkrFitTrack: public TkrRecInfo, virtual public ContainedObject
{    
public:
    /// Constructor/destructor for the class
    TkrFitTrack();
   ~TkrFitTrack();

    void     initializeInfo(unsigned int xgaps, unsigned int ygaps, 
        unsigned int x1st, unsigned int y1st);
    void     initializeQual(double chiSq, double ChiSqSmooth, 
        double rms, double quality, double e, double ms);



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
    
    /// Utilities 
    void   clear();
    void   writeOut(MsgStream& log) const; 

    /// Access to primary quantities on track quality and scattering info
    inline double getChiSquare()           const {return m_chisq;}
    inline double getChiSquareSmooth()     const {return m_chisqSmooth;}
    inline double getScatter()             const {return m_rmsResid;}
    inline double getKalThetaMS()          const {return m_KalThetaMS;}
    inline double getKalEnergy()           const {return m_KalEnergy;}

    /// Access to derived information on gaps
    int           getNumGaps()             const {return m_Xgaps + m_Ygaps;}
    inline int    getNumXGaps()            const {return m_Xgaps;}
    inline int    getNumYGaps()            const {return m_Ygaps;}
    inline int    getNumXFirstGaps ()      const {return m_XistGaps;}
    inline int    getNumYFirstGaps ()      const {return m_YistGaps;}

    /// Access to plane information
    int           getNumHits()             const {return m_hits.size();}
    TkrFitPlane   getFoLPlane(TrackEnd end = Start) const;

    TkrFitPlaneConPtr getHitIterBegin()    const {return m_hits.begin();}
    TkrFitPlaneConPtr getHitIterEnd()      const {return m_hits.end();}

    /// Add hits to our track 
    void          addPlane(TkrFitPlane& newPlane) {m_hits.push_back(newPlane);}

    /// Members below are protected - this gives access to the
    /// classes which fit the tracks (and inherit from this class)
protected:	
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

    // Vector containing hits for each plane traversed
    TkrFitPlaneCol m_hits;
};

//typedef for the Container
  typedef ObjectVector<TkrFitTrack>      TkrFitTrackCol;
  typedef TkrFitTrackCol::const_iterator TkrFitConPtr;
  typedef TkrFitTrackCol::iterator       TkrFitColPtr;

}; //Namespace

#endif
