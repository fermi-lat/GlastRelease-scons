
#include "ParticleTransporter.h"
#include "G4Geantino.hh"
#include "G4StepStatus.hh"
#include "G4VSensitiveDetector.hh"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <string>

//Constructor for the propagator class
ParticleTransporter::ParticleTransporter() 
{
    clearStepInfo();

    return;
}

ParticleTransporter::~ParticleTransporter()
{
    clearStepInfo();
    
    return;
}

void ParticleTransporter::setInitStep(const Point& start, const Vector& dir, const double step)
{
    // Purpose and Method:  Initializes to the starting point and direction and takes initial step
    // Inputs:  The starting point and direction and the initial step length to take
    // Outputs:  None
    // Dependencies: None
    // Restrictions and Caveats: None

    //Clean out any previous step information
    clearStepInfo();

    //Initialize our starting condtions
    startPoint = start;
    startDir   = dir;
    maxArcLen  = step;

    //If step >= 0 then we take the initial step
    if (step >= 0.) transport();

    //Set to negative number so transporter will step to next plane
    maxArcLen  = -1.;

    return;
}

//This method does the track either to the next sensitive plane or to some arc length
bool ParticleTransporter::transport()
{
    // Purpose and Method:  Uses Geant4 to do the transport of the particle from initial to final points
    // Inputs:  None, requires initialization via SetInitStep
    // Outputs:  bool, returns true if target stop point reached, false otherwise
    // Dependencies: None
    // Restrictions and Caveats: None

    bool                  success     = false;
    double                time0       = 0.;
    double                energyPart  = 30.;
    G4ThreeVector         initPoint   = startPoint;

    //If we have already taken a step, start at the last point
    if (stepInfo.size() > 0) initPoint = stepInfo.back()->GetCoords();

    //Set up the Geant tracking classes (define them directly on the stack)
    G4ParticleDefinition* particle    = G4Geantino::Geantino();
    G4DynamicParticle     particleDef(particle, startDir, energyPart);
    G4Track               trackClass(&particleDef, time0, initPoint);
    G4SteppingManager     stepManagerClass;

    G4Track*              track       = &trackClass;
    G4SteppingManager*    stepManager = &stepManagerClass;

    //Set the initial step
    stepManager->SetInitialStep(track);

    //Give the stepping manager the list of processes - for a Geantino there are none
    stepManager->GetProcessNumber();

    //Give the track the pointer to the step
    track->SetStep(stepManager->GetStep());

    //Inform beginning of tracking to physics processes
    track->GetDefinition()->GetProcessManager()->StartTracking();

    //Determine the minimum step to take
    double minStep = minStepSize(stepManager);

    G4StepStatus stepStatus;
    G4double     arcLen      = 0.;
    int          nBoundaries = 0;

    //Record our starting point
    stepInfo.push_back(new TransportStepInfo(track->GetPosition(), 0., stepManager->GetfCurrentVolume()));

    //Continue stepping until the track is no longer "alive"
    while(track->GetTrackStatus() == fAlive)
    {
        track->IncrementCurrentStepNumber();
        stepStatus = stepManager->Stepping();

        //Ok, for a Geantino this should be the only step possibility
        //while inside the "world"
        if (stepStatus == fGeomBoundary)
        {
            double trackLen = track->GetTrackLength();
        
            G4VSensitiveDetector* pSensitive = stepManager->GetfSensitive();

            //Set the stopping condition
            //If stepping to next sensitive boundary, then pSensitive in nonzero
            //If stepping a finite distance, total arclen exceeds that distance
            if ((pSensitive && trackLen > minStep && maxArcLen < 0.) || (trackLen >= maxArcLen && maxArcLen >= 0.)) 
            {
                track->SetTrackStatus(fSuspend);
                success = true;
            }

            nBoundaries++;

            stepInfo.push_back(new TransportStepInfo(track->GetPosition(), trackLen-arcLen, stepManager->GetfCurrentVolume()));

            arcLen = trackLen;
        }
        else track->SetTrackStatus(fSuspend);
    }

    //Let the physics processes have a rest
    track->GetDefinition()->GetProcessManager()->EndTracking();

    //If stepping an arc length, make sure we step the desired distance
    if (maxArcLen > 0 && arcLen > maxArcLen)
    {
        G4ThreeVector hitPos  = stepInfo.back()->GetCoords();
        G4double      stepLen = stepInfo.back()->GetArcLen();
        G4double      corLen  = maxArcLen - arcLen;

        hitPos -= corLen * startDir;

        stepInfo.back()->SetCoords(hitPos);
        stepInfo.back()->SetArcLen(stepLen-corLen);
    }

    //If ending in a sensitive layer then stop in the middle of that volume
    G4VPhysicalVolume* pCurVolume = stepInfo.back()->GetVolume();
    if (maxArcLen < 0 && pCurVolume->GetLogicalVolume()->GetSensitiveDetector())
    {
        G4ThreeVector hitPos = stepInfo.back()->GetCoords();
        G4double      corLen = stepInfo.back()->GetArcLen() * 0.5;

        hitPos -= corLen * startDir;

        stepInfo.back()->SetCoords(hitPos);
        stepInfo.back()->SetArcLen(corLen);
    }

    //Stepping complete
    return success;
}

double ParticleTransporter::minStepSize(G4SteppingManager* stepManager)
{
    // Purpose and Method:  Determine a minimum step size to take. If starting in a
    //                      sensitive volume and wanting to stop at the next sensitive
    //                      volume, then we need to declare a minimum step size which 
    //                      will take us out of the current volume before stopping
    // Inputs:  A pointer to the Geant4 stepping manager
    // Outputs:  a double giving the minimum step length to take
    // Dependencies: None
    // Restrictions and Caveats: None

    double minStep = 0.;

    //If we are starting out in a sensitive volume then the presumption
    //is that we want to step past the edge of this volume to the next volume.
    G4VPhysicalVolume* pCurVolume = stepManager->GetfCurrentVolume();
    if (pCurVolume->GetLogicalVolume()->GetSensitiveDetector())
    {
        const G4VTouchable* touchable = stepManager->GetfTouchable1();

        // this transforms it to local coordinates
        HepTransform3D global(*(touchable->GetRotation()), touchable->GetTranslation());

        HepTransform3D local = global.inverse();

        //G4ThreeVector  startPos = track->GetPosition() - touchable->GetTranslation();
        G4ThreeVector  trackPos = stepManager->GetfTrack()->GetPosition();
        G4ThreeVector  trackDir = stepManager->GetfTrack()->GetMomentumDirection();
        trackPos = local * (HepPoint3D)trackPos;

        //This calculates the distance along the direction of the track to the 
        //boundary of the current volume
        minStep = pCurVolume->GetLogicalVolume()->GetSolid()->DistanceToOut(trackPos,trackDir);
    }

    return minStep;
}


void ParticleTransporter::clearStepInfo()
{
    // Purpose and Method:  Clears the vector of TransportStepInfo objects. Typically before 
    //                      initiating a new step
    // Inputs:  None
    // Outputs:  None
    // Dependencies: None
    // Restrictions and Caveats: None

    StepPtr stepPtr = stepInfo.begin();

    while(stepPtr < stepInfo.end()) delete *stepPtr++;

    stepInfo.clear();

    return;
}
