#ifndef IPropagator_H
#define IPropagator_H

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/TkrRecon/TkrTrackParams.h"
#include "idents/VolumeIdentifier.h"
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include <iostream>

static const InterfaceID IID_IPropagator("IPropagator", 1 , 0);

/** @class IPropagator
* @brief Definition of the interface to the particle propagator
*
*
* @author Tracy Usher
* $Header$
*/
class  IPropagator
{
public:
  	/// Retrieve interface ID
	static const InterfaceID& interfaceID() { return IID_IPropagator; }

    //! Methods to intialize propagator before stepping 
    //! Tracking from intitial point and direction
    virtual void setStepStart(const Point& startPos, const Vector& startDir) = 0;
    //! Tracking from initial parameters with up vs down specified
    virtual void setStepStart(const Event::TkrTrackParams& trackPar, double z, bool upwards = false) = 0;

    //! Takes a step of distance given by arcLen
    virtual void step(double arcLen) = 0;

    //! Number of steps taken
    virtual int getNumberSteps() const = 0;

    //! Methods to return information after the step
    //! Return position after stepping, arcLen can be less than step taken
    virtual Point getPosition(double arcLen = -1.) const = 0;
    //! returns the arc length for this step
    virtual double getArcLen() const = 0;

    //! return position at the center of the step in the volume
    virtual Point getStepPosition(int stepIdx = -1) const = 0;
    //! returns the arc length for this step
    virtual double getStepArcLen(int stepIdx = -1) const = 0;

    //! Return track parameters after stepping, arcLen can be less than step taken
    virtual Event::TkrTrackParams getTrackParams(double arcLen   = -1.,
                                                 double momentum = 1.,
                                                 bool   forward  = true)    const = 0;

    //! Return multiple scattering matrix after stepping, arcLen can be less than step taken
    virtual CLHEP::HepMatrix getMscatCov(double arcLen   = -1.,
                                         double momentum = 1.,
                                         bool   forward  = true)   const = 0;

    //! Return volume identifer after stepping
    virtual idents::VolumeIdentifier getVolumeId(double arcLen = -1.) const = 0;
    //! Return volume identifer at given step index
    virtual idents::VolumeIdentifier getStepVolumeId(int stepIdx = -1) const = 0;

    //! Return material name 
    virtual std::string getMaterialName(int stepIdx = -1)  const = 0;
    //! Return the index into the material properties table
    virtual int         getMaterialIndex(int stepIdx = -1) const = 0;

    //! Return radiation lengths traversed
    virtual double getRadLength(double arcLen = -1.) const = 0;
    //! Return radiation lengths for this step
    virtual double getStepRadLength(int stepIdx = -1.) const = 0;

    //! Return the number of sensitive planes crossed to end of step
    virtual int getNumSensePlanesCrossed() const = 0;

    //! Is the point at the end of stepping inside an Active Area? 
    virtual double isInsideActArea()   const = 0;
    virtual double isInsideActLocalX() const = 0;
    virtual double isInsideActLocalY() const = 0;

    //! Is the plane at the end of stepping an X or a Y plane?
    virtual bool   isXPlane()          const = 0;
    virtual bool   isYPlane()          const = 0;

    //! Is the strip at the end of stepping live?
    virtual bool   isStripLive()       const = 0; 

    /// dump current status, to the stream
    virtual void printOn(std::ostream& str=std::cout )const=0;

    virtual ~IPropagator(){};
protected:
    IPropagator(){};

};

#endif

