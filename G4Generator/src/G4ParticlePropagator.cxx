
#include "G4ParticlePropagator.h"
#include "RunManager.h"
#include "DetectorConstruction.h"
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <string>

// static pointer to the singleton object
G4ParticlePropagator* G4ParticlePropagator::m_instance = 0;

G4ParticlePropagator* G4ParticlePropagator::instance()
{
    // Purpose and Method:  returns a pointer to the object, instantiates if it doesn't exist
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

    //We need to retrieve a pointer to the idents map
    RunManager* pRunManager = RunManager::GetRunManager();

    //If the run manager has been instantiated, then we get the map
    if (pRunManager)
    {
        //Use the run manager to look up the detector volume-idents map
        const G4VUserDetectorConstruction* pUser     = pRunManager->GetUserDetectorConstruction();
        const DetectorConstruction*        pDet      = dynamic_cast<const DetectorConstruction*>(pUser);
        DetectorConstruction*              pDetector = const_cast<DetectorConstruction*>(pDet);

        //Finally...
        m_IdMap = pDetector->idMap();
    }
    else
    //Else we have to figure some failure condition
    {
        m_IdMap = 0;
    }

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

//Drives the tracking to go to the next sensitive plane of the same type as the start
//**Not currently implemented, currently does the same as trackToNextPlane()
bool G4ParticlePropagator::trackToNextSamePlane()
{
    // Purpose and Method:  Interface to actual track propagation method
    // Inputs:  None
    // Outputs:  bool, true if successfully reached target plane, false otherwise
    // Dependencies: Must call setStepStart first
    // Restrictions and Caveats: None
    return transport();
}

//Starting the the current volume, this builds the complete volume identifier string 
//needed to specify where we are in the GLAST world
idents::VolumeIdentifier G4ParticlePropagator::constructId(G4VPhysicalVolume* pVolume) const
{
    // Purpose and Method:  Constructs the complete volume indentifier from current volume
    // Inputs:  A pointer to a G4VPhysicalVolume object giving the current (starting) volume
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
    // Purpose and Method:  Returns the number of planes crossed. Tricky, first volume is the
    //                      start point. The first swim is to this boundary, so the volume is
    //                      repeated. So, number of boundaries is size of vector - 2
    // Inputs:  None
    // Outputs:  Integer number of planes crossed
    // Dependencies: None
    // Restrictions and Caveats: None

    return getNumberSteps() > 2 ? getNumberSteps() - 2 : 0;
}

//Return the position at the end of tracking
Point G4ParticlePropagator::position() const
{
    // Purpose and Method:  Returns the position of the track at the end of tracking
    // Inputs:  None
    // Outputs:  a Point representing the final position
    // Dependencies: None
    // Restrictions and Caveats: None

    Point final(0.,0.,0.);

    if (getNumberSteps() > 0)
    {
        G4ThreeVector  stopPoint = getLastStep()->GetCoords();

        final = Point(stopPoint.x(),stopPoint.y(),stopPoint.z());
    }

    return final;
}

//How far did we go?
double G4ParticlePropagator::arcLength() const
{
    // Purpose and Method:  Returns the arc length subtended at the end of tracking
    // Inputs:  None
    // Outputs:  a double representing the total arc length 
    // Dependencies: None
    // Restrictions and Caveats: None

    double totArcLen = 0.;

    ConstStepPtr stepPtr = getStepStart();

    // We keep track of the step length in each volume, so must loop through and add
    while(stepPtr < getStepEnd())
    {
        TransportStepInfo* curStep = *stepPtr++;

        totArcLen += curStep->GetArcLen();
    }

    return totArcLen;
}

//What is the total number of radiation lengths encountered?
double G4ParticlePropagator::radLength() const 
{
    // Purpose and Method:  Returns the arc length subtended at the end of tracking
    // Inputs:  None
    // Outputs:  a double representing the total arc length 
    // Dependencies: None
    // Restrictions and Caveats: None

    double radLen = 0.;

    //Set up iterator for stepping through all the layers
    ConstStepPtr stepPtr = getStepStart();

    while(stepPtr < getStepEnd())
    {
        TransportStepInfo* curStep = *stepPtr++;

        G4VPhysicalVolume* pCurVolume = curStep->GetVolume();
        G4Material*        pMaterial  = pCurVolume->GetLogicalVolume()->GetMaterial();

        double matRadLen = pMaterial->GetRadlen();
        if (matRadLen > 0.) radLen += 1. / matRadLen;
    }

    if (radLen > 0.) radLen = 1. / radLen;

    return radLen;
}

//Are we inside the active area?
//*** NOT YET IMPLEMENTED ***
float G4ParticlePropagator::insideActArea() const
{
    // Purpose and Method:  Returns the distance to the nearest boundary of the
    //                      nearest sensitive volume. 
    // Inputs:  None
    // Outputs:  distance to the edge, >=0 if inside, < 0 if outside
    // Dependencies: None
    // Restrictions and Caveats: None

    G4VPhysicalVolume* pCurVolume = getLastStep()->GetVolume();

    return 0.;
}

// Is the current plane an X plane
// ** This should be replaced with a routine which returns the volume identifier **
bool G4ParticlePropagator::isXPlane() const
{
    // Purpose and Method:  Determines if current stop point is an "X" plane
    // Inputs:  None
    // Outputs:  a bool, true if current point is an "X" plane, false otherwise 
    // Dependencies: None
    // Restrictions and Caveats: None

    G4VPhysicalVolume* pCurVolume = getLastStep()->GetVolume();

    idents::VolumeIdentifier id = constructId(pCurVolume);

    int         trayNum    = id[4];
    int         botTop     = id[6];
    int         view       = id[5];
    int         layer      = trayNum - 1 + botTop;

    //std::string vId_string = id.name();

    return view == 0;
}

// method to sum up the multiple scattering contributions to the track covr. matrix
// over a distance s (s must be less than maxArcLen).
// This is the exact method used in KalParticle located in package Recon
// and authored by Bill Atwood.
HepMatrix G4ParticlePropagator::mScat_Covr(float momentum, float s) const
{    
    // Purpose and Method:  Calculates Q, the 4x4 covariance matrix due to multiple scattering
    // Inputs:  The total momentum of the track and the arc length over which to calculate
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
        TransportStepInfo* curStep = *stepPtr++;

        G4VPhysicalVolume* pCurVolume = curStep->GetVolume();
        G4Material*        pMaterial  = pCurVolume->GetLogicalVolume()->GetMaterial();

        float radLengths = pMaterial->GetRadlen();
        float x0s        = 10000000.;
        float s_dist     = curStep->GetArcLen();
        float s_distp    = s_dist;

        if (radLengths > 0.) x0s = s_dist / radLengths;

        if(dist+s_dist > s) { // pro-rate the last step: s_distp
            if(s_dist > 0) x0s *= (s - dist)/s_dist;
            s_distp = s - dist; 
        }

        if(x0s != 0) {
            float ms_Angle = 14.0*sqrt(x0s)*(1+0.038*log(x0s))/momentum; //MeV
            float ms_Dst  = (s - dist - s_distp)*ms_Angle;  // Disp. over remaining traj
            float ms_sDst = s_distp*ms_Angle/1.7320508;     // Disp. within step
            float ms_Dist = ms_Dst*ms_Dst + ms_sDst*ms_sDst;

            scat_dist  += ms_Dist;
            scat_angle += ms_Angle*ms_Angle;
            scat_covr  += sqrt(ms_Dist)*ms_Angle;		  
        }
        dist += s_dist;
        if(dist >= s ) break;
    }

    Vector startDir = getStartDir();
    double slopeX   = startDir.x()/startDir.z(); 
    double slopeY   = startDir.y()/startDir.z();
    float  cosThsqX = 1./(1.+slopeX*slopeX);
    float  cosThX   = sqrt(cosThsqX);
    float  cosThsqY = 1./(1.+slopeY*slopeY);
    float  cosThY   = sqrt(cosThsqY);
    
    // Create, fill and return the covariance matrix
    HepMatrix cov(4,4,0);
    cov(2,2) = scat_angle/(cosThsqX*cosThsqX);
    cov(1,1) = scat_dist/(cosThX*cosThsqX);
    cov(1,2) = cov(2,1) = scat_covr/cosThsqX;
    cov(4,4) = scat_angle/(cosThsqY*cosThsqY);
    cov(3,3) = scat_dist/(cosThY*cosThsqY);
    cov(3,4) = cov(4,3) = scat_covr/cosThsqY;

    return cov;
}


void G4ParticlePropagator::printOn(std::ostream& str )const
{
/*
  GParticle::printOn(str);
  str << " status: ";
  switch( m_status ) {
     case ALIVE:      str << "alive";    break;
     case LEFT:       str << "left";     break;
     case LOST:       str << "lost";     break;
     case STUCK:      str << "stuck";    break;
	 case DONE:       str << "done";     break;
     default:         str << "unknown";  break;
	}
  str <<'\n';
  str <<" Total Step length: "<<m_arcLength<<"  Proper Time: "<<m_properTime<<'\n';
  str <<" Step Count: "<<m_stepcount<<'\n';

  int num_step = m_stepList.size();
  for(int i=0; i<num_step; i++){
	  m_stepList[i].printOn(str);
  }
*/
  return;
}
