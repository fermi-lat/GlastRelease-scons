// File and Version Information:
// $Header$
//
// Description: Geant4 class for particle transport management
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
G4ParticlePropagator::G4ParticlePropagator(): ParticleTransporter(G4PropagatorTool::TransportationManager)
{
  // Purpose and Method:  Instantiates if it doesn't exist
  // Inputs:  None
  // Outputs:  None
  // Dependencies: Requires that the Geant4 Run Manager has been instantiated
  // Restrictions and Caveats:  See above

  m_IdMap = G4PropagatorTool::VolIdentMap;

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
  // Inputs:  None
  // Outputs:  bool, true if successfully reached target plane, false otherwise
  // Dependencies: Must call setStepStart first
  // Restrictions and Caveats: 

  return transport();
}

//Drives the tracking to go to the next sensitive plane of the same type as the
//start **Not currently implemented, currently does the same as
//trackToNextPlane()
bool G4ParticlePropagator::trackToNextSamePlane()
{
  // Purpose and Method:  Interface to actual track propagation method
  // Inputs:  None
  // Outputs:  bool, true if successfully reached target plane, false otherwise
  // Dependencies: Must call setStepStart first
  // Restrictions and Caveats: None
  return transport();
}

//Starting the the current volume, this builds the complete volume identifier
//string needed to specify where we are in the GLAST world
idents::VolumeIdentifier G4ParticlePropagator::constructId(G4VPhysicalVolume* pVolume) const
{
  // Purpose and Method: Constructs the complete volume indentifier from current
  //                     volume
  // Inputs: A pointer to a G4VPhysicalVolume object giving the current
  //         (starting) volume
  // Outputs:  a VolumeIdentifier
  // Dependencies: Requires that the volume - idents map has been found
  // Restrictions and Caveats: None
  
  using  idents::VolumeIdentifier;
  VolumeIdentifier ret;

  // Loop through volumes until there is no longer a mother volume
  while(pVolume->GetMother())
    {
      // Look up the identifier for this volume
      VolumeIdentifier id = (*m_IdMap)[pVolume];

      // Add this volume's identifier to our total id
      ret.prepend(id);

      // Get next higher volume
      pVolume = pVolume->GetMother();
    }

  return ret;
}

int G4ParticlePropagator::numberPlanesCrossed() const
{
  // Purpose and Method: Returns the number of planes crossed. Tricky, first
  // volume is the //start point. The first swim is to this boundary, so the
  // volume is repeated. So, number of boundaries is size of vector - 2
  // Inputs:  None
  // Outputs:  Integer number of planes crossed
  // Dependencies: None
  // Restrictions and Caveats: None

  return getNumberSteps() > 1 ? getNumberSteps() - 1 : 0;
}

//Return the position at the end of tracking
Point G4ParticlePropagator::position() const
{
  // Purpose and Method: Returns the position of the track at the end of
  //                     tracking
  // Inputs:  None
  // Outputs:  a Point representing the final position
  // Dependencies: None
  // Restrictions and Caveats: None

  Point final(0.,0.,0.);

  if (getNumberSteps() > 0)
    {
      G4ThreeVector  stopPoint = getLastStep().GetCoords();

      final = Point(stopPoint.x(),stopPoint.y(),stopPoint.z());
    }

  return final;
}

//How far did we go?
double G4ParticlePropagator::arcLength() const
{
  // Purpose and Method: Returns the arc length subtended at the end of tracking
  // Inputs:  None
  // Outputs:  a double representing the total arc length 
  // Dependencies: None
  // Restrictions and Caveats: None

  double totArcLen = 0.;

  ConstStepPtr stepPtr = getStepStart();

  // We keep track of the step length in each volume, so must loop through and
  // add
  while(stepPtr < getStepEnd())
    {
//      TransportStepInfo* curStep = *stepPtr++;

//      totArcLen += curStep->GetArcLen();
      totArcLen += (*stepPtr++).GetArcLen();
    }

  return totArcLen;
}

//What is the total number of radiation lengths encountered?
double G4ParticlePropagator::radLength() const 
{
  // Purpose and Method: Returns the arc length subtended at the end of tracking
  // Inputs:  None
  // Outputs:  a double representing the total arc length 
  // Dependencies: None
  // Restrictions and Caveats: None

  return radLength(s);
}

//What is the total number of radiation lengths encountered?
double G4ParticlePropagator::radLength(double arcLen) const 
{
  // Purpose and Method: Returns the arc length subtended at the end of tracking
  // Inputs:  None
  // Outputs:  a double representing the total arc length 
  // Dependencies: None
  // Restrictions and Caveats: None

  double radLen  = 0.;
  double totDist = 0.;

  //Set up iterator for stepping through all the layers
  ConstStepPtr stepPtr = getStepStart();

  while(stepPtr < getStepEnd())
    {
      TransportStepInfo  curStep = *stepPtr++;

      G4VPhysicalVolume* pCurVolume = curStep.GetVolume();
      G4Material*        pMaterial  = pCurVolume->GetLogicalVolume()->GetMaterial();


      double matRadLen = pMaterial->GetRadlen();
      double x0s       = 0.;
      double s_dist    = curStep.GetArcLen();

      if (matRadLen > 0.) x0s = s_dist / matRadLen;

      // Check that next step is within max arc length
      if(totDist+s_dist > arcLen) 
      { 
        // pro-rate the last step: s_distp
        if(s_dist > 0) x0s *= (arcLen - totDist)/s_dist;
        s_dist = arcLen - totDist;
      }

      if (matRadLen > 0.) radLen += x0s;

      if ((totDist += s_dist) > arcLen) break;
    }

  return radLen;
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

//  G4VPhysicalVolume* pCurVolume = getLastStep()->GetVolume();

  double dist = insideActiveArea();

  return dist;
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

  return insideActiveLocalX();
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

  return insideActiveLocalY();
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

  G4VPhysicalVolume* pCurVolume = getLastStep().GetVolume();

  idents::VolumeIdentifier id = constructId(pCurVolume);

  int         trayNum    = id[4];
  int         botTop     = id[6];
  int         view       = id[5];
  int         layer      = trayNum - 1 + botTop;

  //std::string vId_string = id.name();

  return view == 0;
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

  float scat_dist  = 0.;
  float scat_angle = 0.; 
  float scat_covr  = 0.; 
  float dist       = 0.;

  ConstStepPtr stepPtr = getStepStart();

  while(stepPtr < getStepEnd())
    {
      TransportStepInfo  curStep = *stepPtr++;

      G4VPhysicalVolume* pCurVolume = curStep.GetVolume();
      G4Material* pMaterial  = pCurVolume->GetLogicalVolume()->GetMaterial();

      float radLengths = pMaterial->GetRadlen();
      float x0s        = 10000000.;
      float s_dist     = curStep.GetArcLen();
      float s_distp    = s_dist;

      if (radLengths > 0.) x0s = s_dist / radLengths;

      if(dist+s_dist > arcLen) { // pro-rate the last step: s_distp
        if(s_dist > 0) x0s *= (arcLen - dist)/s_dist;
        s_distp = arcLen - dist; 
      }

      if(x0s != 0) {
        float ms_Angle = 14.0*sqrt(x0s)*(1+0.038*log(x0s))/momentum; //MeV
        float ms_Dst  = (arcLen - dist - s_distp)*ms_Angle; // Disp. over remaining traj
        float ms_sDst = s_distp*ms_Angle/1.7320508; // Disp. within step
        float ms_Dist = ms_Dst*ms_Dst + ms_sDst*ms_sDst;

        scat_dist  += ms_Dist;
        scat_angle += ms_Angle*ms_Angle;
        scat_covr  += sqrt(ms_Dist)*ms_Angle;		  
      }
      dist += s_dist;
      if(dist >= arcLen ) break;
    }

  Vector startDir = getStartDir();
  double slopeX = startDir.x()/startDir.z(); 
  double slopeY = startDir.y()/startDir.z();
  double norm_term = 1. + slopeX*slopeX + slopeY*slopeY;

  // The below taken from KalParticle (by Bill Atwood) in order to match results
  // Calc, the matrix elements (see Data Analysis Tech. for HEP, Fruhwirth et al)
  double p33 = (1.+slopeX*slopeX)*norm_term;
  double p34 = slopeX*slopeY*norm_term;
  double p44 = (1.+slopeY*slopeY)*norm_term; 

  //Go from arc-length to Z 
  scat_dist /=  norm_term;
  scat_covr  /=  sqrt(norm_term);

  HepMatrix cov(4,4,0);
  cov(1,1) = scat_dist*p33;
  cov(2,2) = scat_angle*p33; 
  cov(3,3) = scat_dist*p44;
  cov(4,4) = scat_angle*p44;
  cov(1,2) = cov(2,1) = -scat_covr*p33;
  cov(1,3) = cov(3,1) = scat_dist*p34;
  cov(1,4) = cov(2,3) = cov(3,2) = cov(4,1) = -scat_covr*p34;
  cov(2,4) = cov(4,2) = scat_angle*p34;
  cov(3,4) = cov(4,3) = -scat_covr*p44; 
  
  return cov;
}


void G4ParticlePropagator::printOn(std::ostream& str )const
{
  str << "\n";

  printStepInfo(str);

  return;
}
