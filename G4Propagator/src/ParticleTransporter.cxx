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
#include "globals.hh"
#include "g4std/vector"

#include <string>
#include <algorithm>

//Constructor for the propagator class
ParticleTransporter::ParticleTransporter(const G4TransportationManager* TransportationManager) 
{
  m_TransportationManager = TransportationManager;

  clearStepInfo();

  return;
}

ParticleTransporter::~ParticleTransporter()
{
  clearStepInfo();
    
  return;
}

void ParticleTransporter::setInitStep(const Point& start,  const Vector& dir)
{
  // Purpose and Method: Initializes to the starting point and direction. This
  //                     will serve to set the initial volume at the start point
  // Inputs: The starting point and direction 
  // Outputs:  None
  // Dependencies: None
  // Restrictions and Caveats: None

  //Clean out any previous step information
  clearStepInfo();

  //Initialize our starting condtions
  startPoint  = start;
  startDir    = dir;

  //Set the start point at the beginning of our list of volumes
  G4Navigator*       navigator = m_TransportationManager->GetNavigatorForTracking();
  G4VPhysicalVolume* pVolume   = navigator->LocateGlobalPointAndSetup(startPoint, 0, false);

  //Record our starting point
  stepInfo.push_back(TransportStepInfo(startPoint, 0., pVolume));

  return;
}

//This method does the track either to the next sensitive plane or to some arc
//length
bool ParticleTransporter::transport(const double step)
{
  // Purpose and Method: Uses Geant4 to do the transport of the particle from
  //                     initial to final points
  // Inputs:  None, requires initialization via SetInitStep
  // Outputs:  bool, returns true if target stop point reached, false otherwise
  // Dependencies: None
  // Restrictions and Caveats: None

  if (step >= 0.) return StepAnArcLength(step);
  else            return StepToNextPlane();
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

  G4Navigator*  navigator      = m_TransportationManager->GetNavigatorForTracking();
  G4bool        RelativeSearch = false;
  bool          success        = false;
  G4ThreeVector curPoint       = startPoint;
  G4ThreeVector curDir         = startDir;
  G4double      arcLen         = 0.;

  //Set the minimum step
  double fudge       = 0.00001; // 0.01 um
  double minimumStep = minStepSize(startPoint, startDir) + fudge;

  //If we have already taken a step, start at the last point
  if (stepInfo.size() > 0) curPoint = stepInfo.back().GetCoords();

  bool trackStatus = true;

  //Continue stepping until the track is no longer "alive"
  while(trackStatus)
    {
      //Use the G4 navigator to locate the current point, then compute the distance
      //to the current volume boundary
      double maxStep = 1000.;
      double safeStep;

      G4VPhysicalVolume* pCurVolume = navigator->LocateGlobalPointAndSetup(curPoint, 0, true, true);
      double             trackLen   = navigator->ComputeStep(curPoint, curDir, maxStep, safeStep);

      //If we are right on the boundary then jump over and recalculate
      if (trackLen == 0.)
      {
        double        fudge      = 0.01; // 10 um over the edge
        G4ThreeVector fudgePoint = curPoint + fudge*curDir;

        pCurVolume = navigator->LocateGlobalPointAndSetup(fudgePoint, 0, true, true);
        trackLen   = navigator->ComputeStep(fudgePoint, curDir, maxStep, safeStep) + fudge;
      }

      //If infinite trackLen then we have stepped outside of GLAST
      if (trackLen > 10000000.) 
      {
        break;
      }

      //Which volume are we currently in?
      G4VPhysicalVolume* pSiVolume  = findSiLadders(pCurVolume);

      // If we found the SiLadders volume, compute the distance to leave
      if (pSiVolume)
      {
          // this transforms it to local coordinates
          G4ThreeVector trackPos = navigator->ComputeLocalPoint(curPoint);
          G4ThreeVector trackDir = navigator->ComputeLocalAxis(curDir);

          trackLen   = pSiVolume->GetLogicalVolume()->GetSolid()->DistanceToOut(trackPos,trackDir);
          pCurVolume = pSiVolume;
      }

      //Where did we end up?
      G4ThreeVector newPoint = curPoint + trackLen * curDir;

      stepInfo.push_back(TransportStepInfo(newPoint, trackLen, pCurVolume));

      arcLen     += trackLen;
      curPoint    = newPoint;
        
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
      G4ThreeVector hitPos = stepInfo.back().GetCoords();
      G4double      corLen = stepInfo.back().GetArcLen() * 0.5;

      hitPos -= corLen * startDir;

      stepInfo.back().SetCoords(hitPos);
      stepInfo.back().SetArcLen(corLen);
    }

  //Stepping complete
  return success;
}


//This method does the track either to the next sensitive plane or to some arc
//length
bool ParticleTransporter::StepAnArcLength(const double maxArcLen)
{
  // Purpose and Method: Uses Geant4 to do the transport of the particle from
  //                     initial to final points
  // Inputs:  None, requires initialization via SetInitStep
  // Outputs:  bool, returns true if target stop point reached, false otherwise
  // Dependencies: None
  // Restrictions and Caveats: None

  G4Navigator*  navigator   = m_TransportationManager->GetNavigatorForTracking();
  bool          success     = false;
  G4ThreeVector curPoint    = startPoint;
  G4ThreeVector curDir      = startDir;
  G4double      arcLen      = 0.;

  //If we have already taken a step, start at the last point
  if (stepInfo.size() > 0) curPoint = stepInfo.back().GetCoords();

  //If maxArcLen is zero then we are done, otherwise we need to track the particle
  if (maxArcLen > 0.)
  {
    bool          trackStatus = true;
    double        fudge       = 0.0;
    G4ThreeVector fudgePoint  = curPoint;

    //Continue stepping until the track is no longer "alive"
    while(trackStatus)
      {
        double maxStep = 1000.;
        double safeStep;

        //Use the G4 navigator to locate the current point, then compute the distance
        //to the volume boundary
        G4VPhysicalVolume* pCurVolume = navigator->LocateGlobalPointAndSetup(fudgePoint, 0, true, true);
        double             trackLen   = navigator->ComputeStep(fudgePoint, curDir, maxStep, safeStep) + fudge;

        //If we are right on the boundary then jump over and recalculate
        if (trackLen == 0.)
        {
          fudge       = 0.01; // 10 um
          fudgePoint += fudge*curDir;

          pCurVolume  = navigator->LocateGlobalPointAndSetup(fudgePoint, 0, true, true);
          trackLen    = navigator->ComputeStep(fudgePoint, curDir, maxStep, safeStep) + fudge;
        }

        //If infinite trackLen then we have stepped outside of GLAST
        if (trackLen > 10000000.) break;

        //Update arc length, point and volume
        arcLen     += trackLen;
        curPoint   += trackLen * curDir;
        fudgePoint  = curPoint + fudge * curDir;
        
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
        stepInfo.push_back(TransportStepInfo(curPoint, trackLen, pCurVolume));
      }
  }

  //Stepping complete
  return success;
}

G4VPhysicalVolume* ParticleTransporter::findSiLadders(G4VPhysicalVolume* pVolume) const
{
    G4VPhysicalVolume* volume = pVolume;

    while(volume)
    {
        G4String volName = volume->GetName();

        if (volName.contains("SiLadders"))break;

        volume = volume->GetMother();
    }

    return volume;
}

double ParticleTransporter::minStepSize(const Point& start, const Vector& dir) const
{
  // Purpose and Method: This function determines the distance from the current point
  //                     (given by start), in the direction dir, to exit an SiLadders
  //                     volume, assuming that the current volume is contained within
  //                     an SiLadders volume. 
  //                     Active silicon is in a volume called SiWaferActive which are
  //                     ultimately contained within volumes whose name contains the
  //                     substring SiLadders. So, this routine will give the distance,
  //                     along direction dir, to step out of a current active volume.
  // Inputs:  The current point and direction
  // Outputs:  a double giving the distance to exit the current volume in the 
  //           direction given by dir
  // Dependencies: None
  // Restrictions and Caveats: None

  double   minStep = 0.;

  G4Navigator*         navigator  = m_TransportationManager->GetNavigatorForTracking();
  G4VPhysicalVolume*   pCurVolume = navigator->LocateGlobalPointAndSetup(start, 0, true, true);
  G4TouchableHistory*  pTouchHis  = navigator->CreateTouchableHistory();

  const G4NavigationHistory* pNavHis    = pTouchHis->GetHistory();

  const G4ThreeVector testPoint = navigator->GetCurrentLocalCoordinate();

  for (int depth = pNavHis->GetDepth(); depth > 0; depth--)
  {
      pCurVolume = pNavHis->GetVolume(depth);

      if (pCurVolume->GetName().contains("SiLadders"))
      {
          const G4AffineTransform& transform = pNavHis->GetTransform(depth);
          G4ThreeVector            trackPos = transform.TransformPoint(start);
          G4ThreeVector            trackDir = transform.IsRotated() 
                                            ? transform.TransformAxis(dir) : dir ;

          minStep = pCurVolume->GetLogicalVolume()->GetSolid()->DistanceToOut(trackPos,trackDir);

          break;
      }
  }

  // Clean up the temporary history we made
  delete pTouchHis;

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
    G4Navigator*  navigator = m_TransportationManager->GetNavigatorForTracking();
    G4ThreeVector stopPoint = (stepInfo.back()).GetCoords();
    Point         curPoint  = Point(stopPoint.x(),stopPoint.y(),stopPoint.z());
    G4double      arcLen    = 0.;
    Vector        curDir    = dir;

    G4VPhysicalVolume* pCurVolume = navigator->LocateGlobalPointAndSetup(curPoint, 0, true, true);

    //If we are already in a sensitive volume then simple calculation for distance
    if (pCurVolume->GetLogicalVolume()->GetSensitiveDetector())
    {
      // this transforms it to local coordinates
      G4ThreeVector trackPos = navigator->ComputeLocalPoint(curPoint);
      G4ThreeVector trackDir = navigator->ComputeLocalAxis(curDir);

      //This calculates the distance along the direction of the track to the 
      //boundary of the current volume
      distToActArea = pCurVolume->GetLogicalVolume()->GetSolid()->DistanceToOut(trackPos,trackDir);
    }
    //Otherwise, we need to search for the nearest active volume
    else
    {
        G4ThreeVector fudgePoint  = curPoint;
        double        arcLen      = 0.;
        double        maxStep     = minStepSize(curPoint, curDir);
        double        fudge       = 0.;
        bool          trackStatus = maxStep > 0.;

        while(trackStatus)
        {
            double safeStep;

            //Use the G4 navigator to locate the current point, then compute the distance
            //to the current volume boundary
            G4VPhysicalVolume* pCurVolume = navigator->LocateGlobalPointAndSetup(fudgePoint, 0, true, true);
            double             trackLen   = navigator->ComputeStep(fudgePoint, curDir, maxStep, safeStep) + fudge;

            //If we are right on the boundary then jump over and recalculate
            if (trackLen == 0.)
            {
                fudge       = 0.01;           // 10 um
                fudgePoint += fudge * curDir;

                pCurVolume  = navigator->LocateGlobalPointAndSetup(fudgePoint, 0, true, true);
                trackLen    = navigator->ComputeStep(fudgePoint, curDir, maxStep, safeStep) + fudge;
            }

            //If we have reached the max step then we're not finding anything
            if (arcLen+trackLen >= maxStep) 
            {
                arcLen = 100000.;
                break;
            }

            // If we found the SiLadders volume, compute the distance to leave
            if (pCurVolume->GetLogicalVolume()->GetSensitiveDetector())
            {
                break;
            }

            //Update distance stepped
            arcLen     += trackLen;
            curPoint    = curPoint + trackLen * curDir;
            fudgePoint  = curPoint + fudge * curDir;
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

//  while(stepPtr < stepInfo.end()) delete *stepPtr++;

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
      TransportStepInfo  stepInfo   = *iter++;
      G4VPhysicalVolume* pCurVolume = stepInfo.GetVolume();
      G4ThreeVector      position   = stepInfo.GetCoords();
      G4Material*        pMaterial  = pCurVolume->GetLogicalVolume()->GetMaterial();
      double             matRadLen  = pMaterial->GetRadlen();
      double             x0s        = 0.;
      double             stepDist   = stepInfo.GetArcLen();

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