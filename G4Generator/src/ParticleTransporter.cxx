// File and Version Information:
// $Header$
//
// Description: Geant4 class for particle transport management
//
// Author(s):
//      T.Usher

#include "ParticleTransporter.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "G4TransportationManager.hh"

#include <string>
#include <algorithm>

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

void ParticleTransporter::setInitStep(const Point& start, 
                                      const Vector& dir, const double step)
{
  // Purpose and Method: Initializes to the starting point and direction and
  //                     takes initial step
  // Inputs: The starting point and direction and the initial step length to
  //         take
  // Outputs:  None
  // Dependencies: None
  // Restrictions and Caveats: None

  //Clean out any previous step information
  clearStepInfo();

  //Initialize our starting condtions
  startPoint  = start;
  startDir    = dir;
  maxArcLen   = step;

  //Set the start point at the beginning of our list of volumes
  G4Navigator*          navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4TouchableHistory    touchHist;
  
  navigator->LocateGlobalPointAndUpdateTouchable(startPoint, &touchHist, true);

  //Record our starting point
  stepInfo.push_back(new TransportStepInfo(startPoint, 0., touchHist.GetVolume()));

  //If step >= 0 then stepping through arc length (and doing it now)
  if (step > 0.) transport();

  //Set the minimum step
  double fudge = 0.00001; // 0.01 um
  minimumStep = minStepSize(start, dir) + fudge;

  //Leave in a state to step to the next sensitive plane
  maxArcLen   = -1.;

  return;
}

//This method does the track either to the next sensitive plane or to some arc
//length
bool ParticleTransporter::transport()
{
  // Purpose and Method: Uses Geant4 to do the transport of the particle from
  //                     initial to final points
  // Inputs:  None, requires initialization via SetInitStep
  // Outputs:  bool, returns true if target stop point reached, false otherwise
  // Dependencies: None
  // Restrictions and Caveats: None

  if (maxArcLen >= 0.) return StepAnArcLength();
  else                 return StepToNextPlane();
}


//This method does the track either to the next sensitive plane or to some arc
//length
bool ParticleTransporter::StepToNextPlane()
{
  // Purpose and Method: Uses Geant4 to do the transport of the particle from
  //                     initial to final points
  // Inputs:  None, requires initialization via SetInitStep
  // Outputs:  bool, returns true if target stop point reached, false otherwise
  // Dependencies: None
  // Restrictions and Caveats: None

  G4Navigator*  navigator   = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  bool          success     = false;
  G4ThreeVector curPoint    = startPoint;
  G4ThreeVector curDir      = startDir;
  G4double      arcLen      = 0.;

  G4VPhysicalVolume* pCurVolume;
  G4TouchableHistory newTouch;

  //If we have already taken a step, start at the last point
  if (stepInfo.size() > 0) curPoint = stepInfo.back()->GetCoords();

  //If maxArcLen is not less than zero then we are done
  if (maxArcLen < 0.)
  {
    bool trackStatus = true;

    //Continue stepping until the track is no longer "alive"
    while(trackStatus)
      {
        //Use the G4 navigator to locate the current point, then compute the distance
        //to the current volume boundary
        double maxStep = 1000.;
        double safeStep;
        navigator->LocateGlobalPointAndUpdateTouchable(curPoint, &newTouch, true);
        double trackLen = navigator->ComputeStep(curPoint, curDir, maxStep, safeStep);

        //If we are right on the boundary then jump over and recalculate
        if (trackLen == 0.)
        {
          double        fudge      = 0.01; // 10 um over the edge
          G4ThreeVector fudgePoint = curPoint + fudge*curDir;

          navigator->LocateGlobalPointAndUpdateTouchable(fudgePoint, &newTouch, true);
          trackLen = navigator->ComputeStep(fudgePoint, curDir, maxStep, safeStep) + fudge;
        }

        //If infinite trackLen then we have stepped outside of GLAST
        if (trackLen > 10000000.) 
        {
          break;
        }

        //Which volume are we currently in?
        G4VPhysicalVolume* pNewVolume = newTouch.GetVolume();
        G4VPhysicalVolume* pSiVolume  = findSiLadders(&newTouch);

        // If we found the SiLadders volume, compute the distance to leave
        if (pSiVolume)
        {
            // this transforms it to local coordinates
            HepTransform3D global(*(newTouch.GetRotation()), 
                                    newTouch.GetTranslation());

            HepTransform3D local = global.inverse();

            //G4ThreeVector  startPos = track->GetPosition() - 
            //touchable->GetTranslation();
            G4ThreeVector  trackPos = curPoint;
            G4ThreeVector  trackDir = curDir;
            trackPos   = local * (HepPoint3D)trackPos;
            trackDir   = local * (HepVector3D)trackDir;
            trackLen   = pSiVolume->GetLogicalVolume()->GetSolid()->DistanceToOut(trackPos,trackDir);
            pNewVolume = pSiVolume;
        }

        //Where did we end up?
        G4ThreeVector newPoint = curPoint + trackLen * curDir;

        stepInfo.push_back(new TransportStepInfo(newPoint, trackLen, pNewVolume));

        arcLen     += trackLen;
        curPoint    = newPoint;
        pCurVolume  = pNewVolume;
        
        //Stop if we have found a new layer
        if (arcLen > minimumStep && pCurVolume->GetName().contains("SiLadders"))
        {
            success = true;
            break;
        }
      }

    //If ending in a sensitive layer then stop in the middle of that volume
    if (success)
      {
        G4ThreeVector hitPos = stepInfo.back()->GetCoords();
        G4double      corLen = stepInfo.back()->GetArcLen() * 0.5;

        hitPos -= corLen * startDir;

        stepInfo.back()->SetCoords(hitPos);
        stepInfo.back()->SetArcLen(corLen);
      }
  }

  //Stepping complete
  return success;
}


//This method does the track either to the next sensitive plane or to some arc
//length
bool ParticleTransporter::StepAnArcLength()
{
  // Purpose and Method: Uses Geant4 to do the transport of the particle from
  //                     initial to final points
  // Inputs:  None, requires initialization via SetInitStep
  // Outputs:  bool, returns true if target stop point reached, false otherwise
  // Dependencies: None
  // Restrictions and Caveats: None

  G4Navigator*  navigator   = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  bool          success     = false;
  G4ThreeVector curPoint    = startPoint;
  G4ThreeVector curDir      = startDir;
  G4double      arcLen      = 0.;

  G4VPhysicalVolume* pCurVolume;
  G4TouchableHistory newTouch;

  //If we have already taken a step, start at the last point
  if (stepInfo.size() > 0) curPoint = stepInfo.back()->GetCoords();

  //If maxArcLen is zero then we are done, otherwise we need to track the particle
  if (maxArcLen > 0.)
  {
    bool trackStatus = true;

    //Continue stepping until the track is no longer "alive"
    while(trackStatus)
      {
        //Use the G4 navigator to locate the current point, then compute the distance
        //to the volume boundary
        double maxStep = 1000.;
        double safeStep;
        navigator->LocateGlobalPointAndUpdateTouchable(curPoint, &newTouch, true);
        double trackLen = navigator->ComputeStep(curPoint, curDir, maxStep, safeStep);

        //If we are right on the boundary then jump over and recalculate
        if (trackLen == 0.)
        {
          double        fudge      = 0.01; // 10 um over the edge
          G4ThreeVector fudgePoint = curPoint + fudge*curDir;

          navigator->LocateGlobalPointAndUpdateTouchable(fudgePoint, &newTouch, true);
          trackLen = navigator->ComputeStep(fudgePoint, curDir, maxStep, safeStep) + fudge;
        }

        //If infinite trackLen then we have stepped outside of GLAST
        if (trackLen > 10000000.) break;

        //Update arc length, point and volume
        arcLen     += trackLen;
        curPoint   += trackLen * curDir;
        pCurVolume  = newTouch.GetVolume();
        
        //Stop condition is that we have reached or overrun our arclength
        if (arcLen >= maxArcLen)
        {
            success     = true;
            trackStatus = false;

            // Make sure we didn't go too far
            if (arcLen > maxArcLen)
            {
                G4double corLen = arcLen - maxArcLen;

                curPoint -= corLen * curDir;
                trackLen -= corLen;
            }
        }

        //Store current step info
        stepInfo.push_back(new TransportStepInfo(curPoint, trackLen, pCurVolume));
      }
  }

  //Stepping complete
  return success;
}

G4VPhysicalVolume* ParticleTransporter::findSiLadders(const G4TouchableHistory* pHistory) const
{
    G4VPhysicalVolume* volume     = 0;
    int                touchDepth = pHistory->GetHistoryDepth();


    //Ok, ugliness to actually be stepping in an SiLadders volume
    //Only going to be in something like this if we have a deep enough history
    if (touchDepth > 5)
    {
        G4String SiLadders("SiLadders");

        // Search for an enclosing SiLadders volume, could be 4 up from active wafer
        for(int idx = 0; idx < 5; idx++)
        {
            G4VPhysicalVolume* testVol = pHistory->GetVolume(idx);
            G4String           volName = testVol->GetName();

            if (volName.contains("SiLadders"))
            {
                volume = testVol;
                break;
            }
        }
    }

    return volume;
}

double ParticleTransporter::minStepSize(const Point& start, const Vector& dir) const
{
  // Purpose and Method: Determine a minimum step size to take. If starting in a
  // sensitive volume and wanting to stop at the next sensitive volume, then we
  // need to declare a minimum step size which will take us out of the current
  // volume before stopping
  // Inputs:  A pointer to the Geant4 stepping manager
  // Outputs:  a double giving the minimum step length to take
  // Dependencies: None
  // Restrictions and Caveats: None

  G4String waferNest("SiWaferNest");
  double   minStep = 0.;

  G4Navigator*  navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4TouchableHistory touchable;
  
  navigator->LocateGlobalPointAndUpdateTouchable(start, &touchable, true);
  G4VPhysicalVolume*    pCurVolume = findSiLadders(&touchable);

  //If we are already in an "SiLadders" layer then calculate the minimum step to
  //get out of this volume and search for the next
  if (pCurVolume)
    {
      // this transforms it to local coordinates
      HepTransform3D global(*(touchable.GetRotation()), 
                              touchable.GetTranslation());

      HepTransform3D local = global.inverse();

      //G4ThreeVector  startPos = track->GetPosition() - 
      //touchable->GetTranslation();
      G4ThreeVector  trackPos = start;
      G4ThreeVector  trackDir = dir;
      trackPos = local * (HepPoint3D)trackPos;
      trackDir = local * (HepVector3D)trackDir;

      //This calculates the distance along the direction of the track to the 
      //boundary of the current volume
      minStep  = pCurVolume->GetLogicalVolume()->GetSolid()->DistanceToOut(trackPos,trackDir);
  }

  return minStep;
}

//Method to return distance to nearest active edge in local X direction
double ParticleTransporter::insideActiveLocalX() const
{
  // Purpose and Method: Finds the distance to the edge of 
  //           of the nearest active area in the X (measurement) direction.
  //           Computes the distance to the nearest edge in pos/neg X direction 
  //           and then returns the shortest distance.
  // Inputs:  None
  // Outputs:  double, the distance to the edge of the nearest active area. If the
  //           last tracked point is inside an active area in that coordinate, 
  //           then the value returned is positive, otherwise it is negative
  // Dependencies: Must have been set up with at least a call to setInitStep.
  // Restrictions and Caveats: None

    Vector unitX(1.,0.,0.);
    return distanceToClosestEdge(unitX);
}

//Finds the distance to the closest active edge in a given direction and its opposite
double ParticleTransporter::distanceToClosestEdge(const Vector& dir) const
{
  // Purpose and Method: Finds the distance to the edge of 
  //           of the nearest active area in both directions along dir and
  //           returns the shortest distance.
  // Inputs:  Direction along which to find distance
  // Outputs:  double, the distance to the edge of the nearest active area. If the
  //           last tracked point is inside an active area in that coordinate, 
  //           then the value returned is positive, otherwise it is negative
  // Dependencies: Must have been set up with at least a call to setInitStep.
  // Restrictions and Caveats: None

    
    double distInPos = distanceToEdge(dir);
    double distInNeg = distanceToEdge(-dir);
    return fabs(distInPos) < fabs(distInNeg) ? distInPos : distInNeg;
}

//Method to return distance to nearest active edge in local Y direction
double ParticleTransporter::insideActiveLocalY() const
{
  // Purpose and Method: Finds the distance to the edge of 
  //           of the nearest active area in the Y (measurement) direction.
  //           Computes the distance to the nearest edge in pos/neg Y direction 
  //           and then returns the shortest distance.
  // Inputs:  None
  // Outputs:  double, the distance to the edge of the nearest active area. If the
  //           last tracked point is inside an active area in that coordinate, 
  //           then the value returned is positive, otherwise it is negative
  // Dependencies: Must have been set up with at least a call to setInitStep.
  // Restrictions and Caveats: None

    Vector unitY(0.,1.,0.);
    return distanceToClosestEdge(unitY);
}


//Method to drive the determination of the distance to the edge of the 
//nearest active area. 
double ParticleTransporter::insideActiveArea() const
{
  // Purpose and Method: Drives the calculation of the distance to the edge of 
  //                     of the nearest active area. Computes the distance to 
  //                     the nearest edge in pos/neg x and y directions and then
  //                     returns the shortest distance. For a point outside the 
  //                     active area in both X and Y, the (negative) distance to the nearest
  //                     corner is returned.
  // Inputs:  None
  // Outputs:  double, the distance to the edge of the nearest active area. If the
  //           last tracked point is inside an active area, then the value returned
  //           is positive, otherwise it is negative
  // Dependencies: Must have been set up with at least a call to setInitStep.
  // Restrictions and Caveats: None


    double distInX    = insideActiveLocalX();
    double distInY    = insideActiveLocalY();

    if(distInX>0 || distInY>0) {return std::min(distInX, distInY);}
    else                       {return -sqrt(distInX*distInX + distInY*distInY);}

}


//Method to calculate the distance to the edge of the nearest active area along a
//unit vector given by dir
double ParticleTransporter::distanceToEdge(const Vector& dir) const
{
  // Purpose and Method: Calculates the distance along a vector dir to the edge 
  //                     of the nearest active area. 
  // Inputs:  dir, a vector giving the direction to search for the nearest edge
  // Outputs:  double, the distance to the edge, if the last tracked point is 
  //           inside an active area then the value returned is positive, otherwise
  //           the value returned is negative.
  // Dependencies: Must have been setup with (at least) a call to setInitStep
  // Restrictions and Caveats: None


  double distToActArea = -100000.;

  // Ok, protect against getting called with no points
  if (stepInfo.size() > 0)
  {
    G4Navigator*  navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    G4ThreeVector stopPoint = stepInfo.back()->GetCoords();
    Point         curPoint  = Point(stopPoint.x(),stopPoint.y(),stopPoint.z());
    G4double      arcLen    = 0.;
    Vector        curDir    = dir;

    G4VPhysicalVolume* pCurVolume;
    G4TouchableHistory newTouch;

    navigator->LocateGlobalPointAndUpdateTouchable(curPoint, &newTouch, true);

    pCurVolume = newTouch.GetVolume();

    //If we are already in a sensitive volume then simple calculation for distance
    if (pCurVolume->GetLogicalVolume()->GetSensitiveDetector())
    {
      // this transforms it to local coordinates
      HepTransform3D global(*(newTouch.GetRotation()), 
                              newTouch.GetTranslation());

      HepTransform3D local = global.inverse();

      G4ThreeVector  trackPos = stopPoint;
      G4ThreeVector  trackDir = curDir;
      trackPos = local * (HepPoint3D)trackPos;
      trackDir = local * (HepVector3D)trackDir;

      //This calculates the distance along the direction of the track to the 
      //boundary of the current volume
      distToActArea = pCurVolume->GetLogicalVolume()->GetSolid()->DistanceToOut(trackPos,trackDir);
    }
    //Otherwise, we need to search for the nearest active volume
    else
    {
        double arcLen  = 0.;
        double maxStep = minStepSize(curPoint, curDir);
        double safeStep;
        bool   trackStatus = true;

        while(trackStatus)
        {
            //Use the G4 navigator to locate the current point, then compute the distance
            //to the current volume boundary
            navigator->LocateGlobalPointAndUpdateTouchable(curPoint, &newTouch, true);
            double trackLen = navigator->ComputeStep(curPoint, curDir, maxStep, safeStep);

            //If we are right on the boundary then jump over and recalculate
            if (trackLen == 0.)
            {
                double        fudge      = 0.01; // 10 um over the edge
                G4ThreeVector fudgePoint = curPoint + fudge*curDir;

                navigator->LocateGlobalPointAndUpdateTouchable(fudgePoint, &newTouch, true);
                trackLen = navigator->ComputeStep(fudgePoint, curDir, maxStep, safeStep) + fudge;
            }

            //If we have reached the max step then we're not finding anything
            if (arcLen+trackLen >= maxStep) 
            {
                arcLen = 100000.;
                break;
            }

            //Which volume are we currently in?
            G4VPhysicalVolume* pNewVolume = newTouch.GetVolume();

            // If we found the SiLadders volume, compute the distance to leave
            if (pNewVolume->GetLogicalVolume()->GetSensitiveDetector())
            {
                break;
            }

            //Update distance stepped
            arcLen   += trackLen;
            curPoint  = curPoint + trackLen * curDir;
        }

        distToActArea = -arcLen;
    }

  }

  return distToActArea;
}


void ParticleTransporter::clearStepInfo()
{
  // Purpose and Method: Clears the vector of TransportStepInfo
  // objects. Typically before initiating a new step
  // Inputs:  None
  // Outputs:  None
  // Dependencies: None
  // Restrictions and Caveats: None

  StepPtr stepPtr = stepInfo.begin();

  while(stepPtr < stepInfo.end()) delete *stepPtr++;

  stepInfo.clear();

  return;
}

 
void ParticleTransporter::printStepInfo(std::ostream& str) const
{
  // Purpose and Method: Provides method for formatted outputting of tracking results
  // Inputs:  a reference to a standard output stream
  // Outputs:  None
  // Dependencies: None
  // Restrictions and Caveats: None
  str << "** Particle at: " << startPoint << ", dir: " << startDir << "\n";
  str << "   Step Count:  "  <<getNumberSteps() <<'\n';

  ConstStepPtr iter    = getStepStart();
  double       arcLen  = 0.;
  int          stepCnt = 0;

  while(iter < getStepEnd())
  {
      TransportStepInfo* stepInfo   = *iter++;
      G4VPhysicalVolume* pCurVolume = stepInfo->GetVolume();
      G4ThreeVector      position   = stepInfo->GetCoords();
      G4Material*        pMaterial  = pCurVolume->GetLogicalVolume()->GetMaterial();
      double             matRadLen  = pMaterial->GetRadlen();
      double             x0s        = 0.;
      double             stepDist   = stepInfo->GetArcLen();

      if (matRadLen > 0.) x0s = stepDist / matRadLen;

      arcLen += stepDist;

      G4String           volName = printVolName(pCurVolume);

      str << "   Step " << ++stepCnt << ", Volume: " << volName << "\n" 
          << ", position: " << position << ", step: " << stepDist << ", radlens: " << x0s << "\n";
                
  }
  str <<"   Total Step length: "<< arcLen <<"\n";

  return;
}

G4String ParticleTransporter::printVolName(const G4VPhysicalVolume* pCurVolume) const
{
    if (pCurVolume)
    {
        return "(" + pCurVolume->GetName() + printVolName(pCurVolume->GetMother()) + ")";
    }
    else return " ";
}