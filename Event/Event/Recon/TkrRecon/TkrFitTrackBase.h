#ifndef __TKRFITTRACKBASE_H
#define __TKRFITTRACKBASE_H 1

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "Event/Recon/TkrRecon/TkrFitPlane.h"
#include "geometry/Ray.h"

/** 
* @class TkrFitTrackBase
*
* @brief TkrRecon Abstract Interface for obtaining "standard" recon information
*
* This class defines the basic TkrRecon track fit output information which 
* the external user is likely to need to perform further analysis. 
* Basically, at one or the other end of the track this includes:
* 1) The track parameters
* 2) The covariance matrix for these parameters
* 3) The z coordinate associated with the above parameters
* 4) The energy estimate of the track
* 5) The layer number the track
* 6) The tower number the track
*
* The concept for this base class was originally proposed by Bill Atwood. 
* The Tracking Software Group has extended the concept into this abstract
* interface as a mechanism for providing the external user (i.e. outside
* TkrRecon) with a single interface point for important information.
*
* @author(s) The Tracking Software Group
*
* $Header$
*/
namespace Event { // Namespace

class TkrFitTrackBase : public TkrFitPlaneCol, virtual public ContainedObject
{
public:
    /// Retrieve track information at the Start or End of the track
    enum TrackEnd { Null, Start, End };

    /// The track quality
    virtual double       getQuality()                       const = 0;

    /// Energy of the track
    virtual double       getEnergy(TrackEnd end = Start)    const = 0;

    /// Layer at this end of the track
    virtual int          getLayer(TrackEnd end = Start)     const = 0;

    /// Tower at this point on the track
    virtual int          getTower(TrackEnd end = Start)     const = 0;
    
    /// Return position at this end of the track
    virtual Point        getPosition(TrackEnd end = Start)  const = 0;
    
    /// Return the direction at this end of the track
    virtual Vector       getDirection(TrackEnd end = Start) const = 0;

    /// Provide a Ray at this end
    virtual Ray          getRay(TrackEnd end = Start)       const = 0;

    /// Return the track parameters at this end of the track
    virtual TkrFitPar    getTrackPar(TrackEnd end = Start)  const = 0;

    /// Return the Z coordinate of the track parameters
    virtual double       getTrackParZ(TrackEnd end = Start) const = 0;

    /// Return the track covariance matrix at this end of the track
    virtual TkrFitMatrix getTrackCov(TrackEnd end = Start)  const = 0;

    /// Utility function to tell us a valid track/vertex exists
    virtual bool         empty(int numHits)                 const = 0;

    /// seems to need a virtual destructor
    virtual ~TkrFitTrackBase() {}
}; 

//typedef for the Container
  typedef ObjectVector<TkrFitTrackBase>      TkrFitTrackCol;
  typedef TkrFitTrackCol::const_iterator     TkrFitConPtr;
  typedef TkrFitTrackCol::iterator           TkrFitColPtr;

}; // Namespace

#endif
