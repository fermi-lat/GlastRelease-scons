#ifndef G4ParticlePropagator_h
#define G4ParticlePropagator_h 
#ifdef WIN32 // for G4 
#include <float.h>
#endif

#include "G4Generator/IG4GeometrySvc.h"
#include "idents/VolumeIdentifier.h"

#include "GlastSvc/Reco/IKalmanParticle.h"
#include "GlastSvc/Reco/IPropagator.h"
#include "ParticleTransporter.h"

#include <vector>

/** 
* @class G4ParticlePropagator
*
* @brief Singleton object providing particle transport (ie error matrix
* extrapolation)
*
* This class implements a particle transport algorithm for extrapolating the
* track error matrix from a start point to and end point by "swimming" through
* the GLAST geometry. It uses the Geant4 instantiation of the geometry and the
* Geant4 tracking routines.  The interface is adapted from that used for
* RCparticle, the gismo propagator.
*
* This version has now been turned into an interface to G4PropagationTool,
* which is interfaced through the IPropagator interface class.
* April 21, 2003 Tracy Usher
*
* @author Tracy Usher
*
*/

//Class definition for the Geant 4 particle propagator
//class G4ParticlePropagator : public IParticlePropagator 
class G4ParticlePropagator : public IKalmanParticle
{
public: 
  /// Method to return pointer to the singleton object (and instantiate if not
  /// already)
    static  G4ParticlePropagator* instance();

    /// Methods for interfacing to the class
    virtual void      setStepStart(const Point& start, 
                                   const Vector& dir, 
                                   const double step) 
    {m_propagator->setStepStart(start,dir); m_propagator->step(step);}
    virtual bool      trackToNextPlane();
    virtual bool      trackToNextSamePlane();
    virtual int       numberPlanesCrossed() const;
    virtual double    insideActArea() const;
    virtual double    insideActLocalX() const;
    virtual double    insideActLocalY() const;
    virtual bool      stripIsLive() const {return true;} //Not implemented
    virtual Point     position()  const             {return m_propagator->getPosition();}
    virtual Point     position(double arcLen) const {return m_propagator->getPosition(arcLen);}
    virtual double    arcLength() const;
    virtual double    radLength() const;
    virtual double    radLength(double arcLen) const;
    virtual bool      isXPlane()  const;
    virtual HepMatrix mScat_Covr(double momentum, double arcLen) const;

    /// dump current status, to the stream
    virtual void printOn(std::ostream& str=std::cout ) const;
  
private:
    /// Private methods
    G4ParticlePropagator();
   ~G4ParticlePropagator();

    /// Private data
    /// Pointer to the class to make it a singleton
    static G4ParticlePropagator* m_instance;

    IPropagator* m_propagator;
};

#endif




