#ifndef G4ParticlePropagator_h
#define G4ParticlePropagator_h 

#include "RunManager.h"
#include "DetectorConstruction.h"
#include "idents/VolumeIdentifier.h"

#include "globals.hh"
#include "g4std/vector"

#include "GlastSvc/Reco/IKalmanParticle.h"
#include "ParticleTransporter.h"

#include <vector>

/** 
* @class G4ParticlePropagator
*
* @brief Singleton object providing particle transport (ie error matrix extrapolation)
*
* This class implements a particle transport algorithm for extrapolating the track 
* error matrix from a start point to and end point by "swimming" through the GLAST
* geometry. It uses the Geant4 instantiation of the geometry and the Geant4 tracking
* routines. 
* The interface is adapted from that used for RCparticle, the gismo propagator.
*
* @author Tracy Usher
*
*/

//Class definition for the Geant 4 particle propagator
//class G4ParticlePropagator : public IParticlePropagator 
class G4ParticlePropagator : public ParticleTransporter, public IKalmanParticle
{
public: 
    /// Method to return pointer to the singleton object (and instantiate if not already)
    static  G4ParticlePropagator* instance();

    /// Methods for interfacing to the class
    virtual void      setStepStart(const Point& start, const Vector& dir, const double step) {setInitStep(start,dir,step);}
    virtual bool      trackToNextPlane();
    virtual bool      trackToNextSamePlane();
    virtual int       numberPlanesCrossed() const;
    virtual float     insideActArea() const;
    virtual bool      stripIsLive() const {return true;} //Not implemented
    virtual Point     position()  const;
    virtual double    arcLength() const;
    virtual double    radLength() const;
    virtual bool      isXPlane()  const;
    virtual HepMatrix mScat_Covr(float momentum, float s) const;

    /// dump current status, to the stream
    virtual void printOn(std::ostream& str=std::cout ) const;
  
private:
    /// Private methods
    G4ParticlePropagator();
   ~G4ParticlePropagator();

    /// Builds the volume identifier starting from the current volume
    idents::VolumeIdentifier constructId(G4VPhysicalVolume* pVolume) const;

    /// Private data
    /// Pointer to the class to make it a singleton
    static G4ParticlePropagator* m_instance;

    /// This is a pointer to the all important volume->idents map
    /// obtained from the RunManager singleton
    DetectorConstruction::IdMap* m_IdMap;
};

#endif




