// File and Version Information:
// $Header$
//
// Description: Geant4 class for particle transport management
//
// April 26, 2003 - Tracy Usher
// This class now acts as an interface to G4PropagationTool. This preserves the old IKalmanParticle
// interface, which is compatible with the gismo propagator. This will eventually be phased out.
//
// Author(s):
//      T.Usher

#include "G4ParticlePropagator.h"
#include "G4PropagatorTool.h"
#include "G4VUserDetectorConstruction.hh"
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <string>

// static pointer to the singleton object
G4ParticlePropagator* G4ParticlePropagator::m_instance = 0;

G4ParticlePropagator* G4ParticlePropagator::instance()
{
  // Purpose and Method: returns a pointer to the object, instantiates if it
  //                     doesn't exist
  // Inputs:  None
  // Outputs:  pointer to the object
  // Dependencies: None
  // Restrictions and Caveats:  None
  if (!m_instance) m_instance = new G4ParticlePropagator();

  return m_instance;
}

//Constructor for the propagator class
G4ParticlePropagator::G4ParticlePropagator()
{
  // Purpose and Method:  Instantiates if it doesn't exist
  // Inputs:  None
  // Outputs:  None
  // Dependencies: Requires that the Geant4 Run Manager has been instantiated
  // Restrictions and Caveats:  See above

  m_propagator = G4PropagatorTool::propagatorTool;

  return;
}

G4ParticlePropagator::~G4ParticlePropagator()
{
  return;
}

//Drives the tracking to go to the next sensitive plane
bool G4ParticlePropagator::trackToNextPlane()
{
  // Purpose and Method:  Interface to actual track propagation method
  // NOTE: This method now does nothing as tracking was performed by call to setup.
  // Inputs:  None
  // Outputs:  bool, true if successfully reached target plane, false otherwise
  // Dependencies: Must call setStepStart first
  // Restrictions and Caveats: 

  return true;
}

//Drives the tracking to go to the next sensitive plane of the same type as the
//start **Not currently implemented, currently does the same as
//trackToNextPlane()
bool G4ParticlePropagator::trackToNextSamePlane()
{
  // Purpose and Method:  Interface to actual track propagation method
  // NOTE: This method now does nothing as tracking was performed by call to setup.
  // Inputs:  None
  // Outputs:  bool, true if successfully reached target plane, false otherwise
  // Dependencies: Must call setStepStart first
  // Restrictions and Caveats: None
  return true;
}

int G4ParticlePropagator::numberPlanesCrossed() const
{
  // Purpose and Method: Returns the number of planes crossed. 
  // Inputs:  None
  // Outputs:  Integer number of planes crossed
  // Dependencies: None
  // Restrictions and Caveats: None

    return m_propagator->getNumSensePlanesCrossed();
}

//How far did we go?
double G4ParticlePropagator::arcLength() const
{
  // Purpose and Method: Returns the arc length subtended at the end of tracking
  // Inputs:  None
  // Outputs:  a double representing the total arc length 
  // Dependencies: None
  // Restrictions and Caveats: None


  return m_propagator->getArcLen() ;
}

//What is the total number of radiation lengths encountered?
double G4ParticlePropagator::radLength() const 
{
  // Purpose and Method: Returns the arc length subtended at the end of tracking
  // Inputs:  None
  // Outputs:  a double representing the total arc length 
  // Dependencies: None
  // Restrictions and Caveats: None

    return m_propagator->getRadLength();
}

//What is the total number of radiation lengths encountered?
double G4ParticlePropagator::radLength(double arcLen) const 
{
  // Purpose and Method: Returns the arc length subtended at the end of tracking
  // Inputs:  None
  // Outputs:  a double representing the total arc length 
  // Dependencies: None
  // Restrictions and Caveats: None

    return m_propagator->getRadLength(arcLen);
}

//Are we inside the active area?
double G4ParticlePropagator::insideActArea() const
{
  // Purpose and Method:  Returns the distance to the nearest boundary of the
  //                      nearest sensitive volume. 
  // Inputs:  None
  // Outputs:  distance to the edge, >=0 if inside, < 0 if outside
  // Dependencies: None
  // Restrictions and Caveats: None

    return m_propagator->isInsideActArea();
}

//Are we inside the active in the X (measurement) direction?
double G4ParticlePropagator::insideActLocalX() const
{
  // Purpose and Method:  Returns the distance to the nearest X boundary of the
  //                      nearest sensitive volume. 
  // Inputs:  None
  // Outputs:  distance to the edge, >=0 if inside, < 0 if outside
  // Dependencies: None
  // Restrictions and Caveats: None

    return m_propagator->isInsideActLocalX();
}

//Are we inside the active in the Y (non-measurement) direction?
double G4ParticlePropagator::insideActLocalY() const
{
  // Purpose and Method:  Returns the distance to the nearest X boundary of the
  //                      nearest sensitive volume. 
  // Inputs:  None
  // Outputs:  distance to the edge, >=0 if inside, < 0 if outside
  // Dependencies: None
  // Restrictions and Caveats: None

    return m_propagator->isInsideActLocalY();
}

// Is the current plane an X plane ** This should be replaced with a routine
// which returns the volume identifier **
bool G4ParticlePropagator::isXPlane() const
{
  // Purpose and Method:  Determines if current stop point is an "X" plane
  // Inputs:  None
  // Outputs:  a bool, true if current point is an "X" plane, false otherwise 
  // Dependencies: None
  // Restrictions and Caveats: None

  return m_propagator->isXPlane();
}

// method to sum up the multiple scattering contributions to the track
// covr. matrix over a distance s (s must be less than maxArcLen).  This is the
// exact method used in KalParticle located in package Recon and authored by
// Bill Atwood.
HepMatrix G4ParticlePropagator::mScat_Covr(double momentum, double arcLen) const
{    
  // Purpose and Method: Calculates Q, the 4x4 covariance matrix due to multiple
  //                     scattering
  // Inputs: The total momentum of the track and the arc length over which to
  //         calculate
  // Outputs:  a 4x4 HepMatrix  
  // Dependencies: To get any answer, must have already stepped
  // Restrictions and Caveats: None

  Event::TkrFitMatrix trackCov = m_propagator->getMscatCov(momentum,arcLen);
    
  HepMatrix cov(4,4,0);
  cov(1,1) = trackCov.getcovX0X0();
  cov(2,2) = trackCov.getcovSxSx(); 
  cov(3,3) = trackCov.getcovY0Y0();
  cov(4,4) = trackCov.getcovSySy();
  cov(1,2) = cov(2,1) = trackCov.getcovX0Sx();
  cov(1,3) = cov(3,1) = trackCov.getcovX0Y0();
  cov(1,4) = cov(2,3) = cov(3,2) = cov(4,1) = trackCov.getcovX0Sy();
  cov(2,4) = cov(4,2) = trackCov.getcovSxSy();
  cov(3,4) = cov(4,3) = trackCov.getcovY0Sy(); 

  return cov;
}


void G4ParticlePropagator::printOn(std::ostream& str )const
{
  m_propagator->printOn(str);

  return;
}
