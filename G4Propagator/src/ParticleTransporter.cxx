// 
// @file ParticleTransporter.cxx
// @brief Source code for the ParticleTransporter Class
// @author Tracy Usher
//
// File and Version Information:
// $Header$
//

#include "ParticleTransporter.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "globals.hh"
#include "geomdefs.hh"
#include "G4StateManager.hh"
#include "G4PropagatorExceptionHandler.h"

#include <stdexcept>
#include <string>
#include <sstream>
#include <algorithm>

#include <iostream>

//Constructor for the propagator class
ParticleTransporter::ParticleTransporter(const G4TransportationManager* TransportationManager,
                                               IG4GeometrySvc* geoSvc) 
{
    // Purpose and Method: Constructor for the ParticleTransporter class
    // Inputs: Pointers to the G4TransportationManager and the G4GeometrySvc interface 
    // Outputs:  None
    // Dependencies: None
    // Restrictions and Caveats: None
    m_TransportationManager = TransportationManager;
    m_geometrySvc           = geoSvc;

    clearStepInfo();

    prevStartPoint  = Point(-10000.,-10000.,-10000.);
    prevStartDir    = Vector(0.,0.,0.);
    prevArcLen      = 0.;
    sameStartParams = false;

    return;
}

ParticleTransporter::~ParticleTransporter()
{
    // Purpose and Method: Destructor for the ParticleTransporter class
    // Inputs: None 
    // Outputs:  None
    // Dependencies: None
    // Restrictions and Caveats: None
    clearStepInfo();
    
    return;
}

void ParticleTransporter::initExceptionHandler()
{
    // Purpose and Method: Set up exception handling for ParticleTransporter class
    // Inputs: none
    // Outputs:  None
    // Dependencies: None
    // Restrictions and Caveats: None

    // Make sure the exception handling is set up for tracking
    G4StateManager*      stateManager = G4StateManager::GetStateManager();
    G4VExceptionHandler* exceptHand   = stateManager->GetExceptionHandler();

    // If no exception handler, register one
    if (exceptHand == 0)
    {
        G4PropagatorExceptionHandler* handler = new G4PropagatorExceptionHandler();
        stateManager->SetExceptionHandler(handler);
    }

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

    // Compare start position and direction to last start and direction
    // This is an attempt to cache info in order to prevent repeated running of propagator
    // over the same parameters
    if (const_cast<Point&>(start) == prevStartPoint && const_cast<Vector&>(dir) == prevStartDir)
    {
        sameStartParams = true;
    }
    // Otherwise, this is a new setup so initialize accordingly
    else
    {
        // Clean out any previous step information
        clearStepInfo();

        // Test for transportation manager!!
        if (!m_TransportationManager) setTransportationManager(m_geometrySvc->getTransportationManager());

        // Set the volume hierarchy for our initial point
        G4Navigator*       navigator = m_TransportationManager->GetNavigatorForTracking();
        G4VPhysicalVolume* pVolume   = navigator->LocateGlobalPointAndSetup(start, 0, true, false);

        // Let's be sure that we are inside a valid volume (ie inside of GLAST - it can happen!)
        if (!pVolume) 
        {
            std::stringstream errorStream;
            errorStream << "ParticleTransporter given invalid initial conditions. pos: " 
                        << start << " dir: " << dir;
            throw std::domain_error(errorStream.str());
        }

        // Get the arclength for the first step, if we can...
        double maxStep  = 1000.;    // Variables used by G4
        double safeStep = 0.01;     // 
        double trackLen = navigator->ComputeStep(start, dir, maxStep, safeStep);

        const G4AffineTransform& globalToLocal = navigator->GetGlobalToLocalTransform();

        //Record our starting point
        stepInfo.push_back(TransportStepInfo(start, dir, trackLen, globalToLocal, pVolume));

        // Update previous state info
        prevStartPoint  = start;
        prevStartDir    = dir;
        prevArcLen      = 0.;
        sameStartParams = false;
    }

    return;
}

//This method does the track either to the next sensitive plane or to some arc
//length
bool ParticleTransporter::transport(const double step)
{
  // Purpose and Method: Uses Geant4 to do the transport of the track from
  //                     initial to final points. 
  // Inputs:  step - arclength to step through. If step is negative, then will
  //          attempt to step until reaching the next sensitive via StepToNextPlane.
  // Outputs:  bool, returns true if target stop point reached, false otherwise
  // Dependencies: Requires initialization via SetInitStep
  // Restrictions and Caveats: None

  if (step >= 0.) return StepAnArcLength(step);
  else            return StepToNextPlane();
}

Point  ParticleTransporter::getStartPoint() const
{
    G4ThreeVector startPoint = getStep(0).GetEntryPoint();
    return Point(startPoint.x(),startPoint.y(),startPoint.z());
}

Vector ParticleTransporter::getStartDir() const
{
    G4ThreeVector startDir = getStep(0).GetDirection();
    return Point(startDir.x(),startDir.y(),startDir.z());
}

//This method does the track either to the next sensitive plane or to some max arc
//length
bool ParticleTransporter::StepToNextPlane()
{
  // Purpose and Method: Uses Geant4 to do the transport of the particle from
  //                     initial to final points
  // ***** This method is no longer supported and remains for historical purposes *****
  // Inputs:  None 
  // Outputs:  bool, returns true if target stop point reached, false otherwise
  // Dependencies: Requires initialization via SetInitStep
  // Restrictions and Caveats: None

  G4Navigator*  navigator = m_TransportationManager->GetNavigatorForTracking();
  bool          success   = false;
  G4ThreeVector curPoint  = stepInfo.back().GetEndPoint();
  G4ThreeVector curDir    = stepInfo.back().GetDirection();
  G4double      arcLen    = stepInfo.back().GetArcLen();

  //Set the minimum step
  double fudge       = 0.00001; // 0.01 um
  double minimumStep = minStepSize(curPoint, curDir) + fudge;

  bool trackStatus = true;

  //Continue stepping until the track is no longer "alive"
  while(trackStatus)
    {
      //Use the G4 navigator to locate the current point, then compute the distance
      //to the current volume boundary
      double maxStep = 1000.;
      double safeStep;

      const G4VPhysicalVolume* pCurVolume    = navigator->LocateGlobalPointAndSetup(curPoint, 0, true, false);
      G4AffineTransform        globalToLocal = navigator->GetGlobalToLocalTransform();
      double                   trackLen      = navigator->ComputeStep(curPoint, curDir, maxStep, safeStep);

      //If we are right on the boundary then jump over and recalculate
      if (trackLen == 0.)
      {
        double        fudge      = 0.01; // 10 um over the edge
        G4ThreeVector fudgePoint = curPoint + fudge*curDir;

        pCurVolume    = navigator->LocateGlobalPointAndSetup(fudgePoint, 0, true, false);
        globalToLocal = navigator->GetGlobalToLocalTransform();
        trackLen      = navigator->ComputeStep(fudgePoint, curDir, maxStep, safeStep) + fudge;
      }

      //If infinite trackLen then we have stepped outside of GLAST
      if (trackLen > 10000000.) 
      {
        break;
      }

      //Which volume are we currently in?
      G4TouchableHistory*      touchable = navigator->CreateTouchableHistory(); 
      const G4VPhysicalVolume* pSiVolume = findSiLadders(touchable);
      delete touchable;

      // If we found the SiLadders volume, compute the distance to leave
      if (pSiVolume)
      {
          // this transforms it to local coordinates
          G4ThreeVector trackPos = globalToLocal.TransformPoint(curPoint);
          G4ThreeVector trackDir = globalToLocal.IsRotated() ? globalToLocal.TransformAxis(curDir) : curDir;

          trackLen   = pSiVolume->GetLogicalVolume()->GetSolid()->DistanceToOut(trackPos,trackDir);
          pCurVolume = pSiVolume;
      }

      //Where did we end up?
      G4ThreeVector newPoint = curPoint + trackLen * curDir;

      stepInfo.push_back(TransportStepInfo(newPoint, curDir, trackLen, globalToLocal, pCurVolume));

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
      G4ThreeVector hitPos = stepInfo.back().GetEndPoint();
      G4double      corLen = stepInfo.back().GetArcLen() * 0.5;

      hitPos -= corLen * curDir;

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
    // Inputs:  None 
    // Outputs:  bool, returns true if target stop point reached, false otherwise
    // Dependencies: Requires initialization via SetInitStep
    // Restrictions and Caveats: None

    // First check if this is a repeated step
    if (sameStartParams && prevArcLen == maxArcLen) return true;

    G4Navigator*  navigator = m_TransportationManager->GetNavigatorForTracking();
    bool          success   = false;
    G4double      arcLen    = stepInfo.back().GetArcLen();
    G4ThreeVector curDir    = stepInfo.back().GetDirection();
    G4ThreeVector curPoint  = stepInfo.back().GetEndPoint();
    int           maxTries  = 3;

    // If arcLen < maxArcLen then we still need to do some tracking
    // Note that this assumes arcLen is the value set at intialization! 
    if (arcLen < maxArcLen)
    {
        bool trackStatus = true;

        //Continue stepping until the track is no longer "alive"
        while(trackStatus)
        {
            double maxStep      = 1000.;    // Variables used by G4
            double safeStep     = 0.1;      // 
            double stepOverDist = 1000. * kCarTolerance;

            // Create a point which we hope is "just over the boundary"
            G4ThreeVector overPoint = curPoint + stepOverDist * curDir;

            //Use the G4 navigator to locate the current point
            G4VPhysicalVolume* pCurVolume    = navigator->LocateGlobalPointAndSetup(overPoint, 0, true, false);

            // Let's be sure that we are inside a valid volume (ie inside of GLAST - it can happen!)
            if (!pCurVolume) 
            {
                std::stringstream errorStream;
                errorStream << "StepAnArcLength cannot find position within GLAST. pos: " 
                    << overPoint << " dir: " << curDir << " step: " << stepInfo.size();
                throw std::domain_error(errorStream.str());
            }

            // If ok, compute the distance to the volume boundary
            G4AffineTransform  globalToLocal = navigator->GetGlobalToLocalTransform();
            double             trackLen      = navigator->ComputeStep(overPoint, curDir, maxStep, safeStep) + stepOverDist;

            //If we are right on the boundary then jump over and recalculate
            while(trackLen <= stepOverDist)
            {
                // We are here because we are "on" the boundary (within the tolerance for that)
                // Two reasons for the problem: 1) travelling almost exactly parallel to the boundary
                // or, 2) not parallel
                // Try to ascertain which here
                G4bool              temp          = false;
                const G4ThreeVector localExitNrml = navigator->GetLocalExitNormal(&temp);
                G4ThreeVector       trackDir      = curDir;

                double trkToExitAng = trackDir.dot(localExitNrml);

                // Parallel track case
                // |cos(track to surface normal)| < 1e-5 - nearly exactly perpendicular
                if (fabs(trkToExitAng) < 10000. * kCarTolerance)
                {
                    // In this case one of the coordinates is either on the surface of the volume or 
                    // within its tolerance (defined by kCarTolerance). So, we try to move the point to
                    // get outside of the tolerance and start over... 

                    overPoint -= kCarTolerance * localExitNrml;
                }
                // Not exactly parallel, here we can probably "help" things by increasing the step by hand
                // and starting again
                else 
                {
                    stepOverDist += 1000. * kCarTolerance / fabs(trkToExitAng); // 1000 * minimum tolerance in G4
                    overPoint    += stepOverDist*curDir;
                }

                // Re-step with new points
                pCurVolume    = navigator->LocateGlobalPointAndSetup(overPoint, 0, true, true);
                globalToLocal = navigator->GetGlobalToLocalTransform();
                trackLen      = navigator->ComputeStep(overPoint, curDir, maxStep, safeStep) + stepOverDist;
            }

            //If infinite trackLen then we have stepped outside of GLAST
            if (trackLen > 10000000.) break;

            //Update arc length
            arcLen += trackLen;
        
            //Stop condition is that we have reached or overrun our arclength
            if (arcLen >= maxArcLen)
            {
                success     = true;
                trackStatus = false;

                // Make sure we didn't go too far
                if (arcLen > maxArcLen)
                {
                    G4double corLen = arcLen - maxArcLen;

                    //curPoint -= corLen * curDir;
                    trackLen -= corLen;
                }
            }

            // Store current information: the point we started the step at, the step length
            // and the pointer to the current G4 logical volume
            stepInfo.push_back(TransportStepInfo(curPoint, curDir, trackLen, globalToLocal, pCurVolume));

            // Now update the current point for the next loop
            curPoint += trackLen * curDir;
        }
    }
    // If the initialization step has overstepped, reset it here
    else stepInfo.back().SetArcLen(maxArcLen);

    // Save step arc length
    prevArcLen = maxArcLen;

    //Stepping complete
    return success;
}

const G4VPhysicalVolume* ParticleTransporter::findSiLadders(G4TouchableHistory* theTouchable) const
{
    // Purpose and Method: Search the volume tree to see if we are inside an "SiLadders" volume
    // Inputs:  Pointer to the starting G4VPhysicalVolume 
    // Outputs:  bool, returns true if target stop point reached, false otherwise
    // Dependencies: Requires initialization via SetInitStep
    // Restrictions and Caveats: None

    const G4VPhysicalVolume* volume = 0;

    if (theTouchable)
	{
	    for( int i = 0; i<theTouchable->GetHistoryDepth() ; ++i) 
        {
	        volume = theTouchable->GetVolume(i); 
            G4String volName = volume->GetName();
            if (volName.contains("SiLadders"))break;
	    }
	}

    return volume;
}

double ParticleTransporter::minStepSize(const G4ThreeVector& start, const G4ThreeVector& dir) const
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

    G4ThreeVector testPoint = navigator->GetGlobalToLocalTransform().TransformPoint(start);
    //**const G4ThreeVector testPoint = navigator->GetCurrentLocalCoordinate();

    for (int depth = pNavHis->GetDepth(); depth > 0; depth--)
    {
        pCurVolume = pNavHis->GetVolume(depth);

        if (pCurVolume->GetName().contains("SiLadders"))
        {
            const G4AffineTransform& transform = pNavHis->GetTransform(depth);
            G4ThreeVector            trackPos = transform.TransformPoint(start);
	        G4ThreeVector trackDir;
	        if (transform.IsRotated()) trackDir = transform.TransformAxis(dir);
	        else                       trackDir = dir;

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

double ParticleTransporter::getTotalArcLen() const
{
    // Purpose and Method: Calculates the total arc length traversed through all steps
    // Inputs:  None
    // Outputs:  Total arc length traversed
    // Dependencies: Must have first tracked a particle through some volumes...
    // Restrictions and Caveats: None

    double arcLen = 0;

    //Set up iterator for stepping through all the layers
    ConstStepPtr stepPtr = getStepStart();

    while(stepPtr < getStepEnd()) arcLen += (*stepPtr++).GetArcLen();

    return arcLen;
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
    else                       {return std::max(distInX, distInY);}
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
        G4ThreeVector curPoint  = stepInfo.back().GetEndPoint();
        G4ThreeVector curDir    = dir;
        double        maxStep   = 200.;

        // This needed to insure the volume hierarchy is setup correctly for steps below
        const G4VPhysicalVolume* pCurVolume = stepInfo.back().GetVolume();
            
        // Transform the global point and direction to local (to this volume)
        const G4AffineTransform& globalToLocal = stepInfo.back().GetGlobalToLocal();
        G4ThreeVector localPos = globalToLocal.TransformPoint(curPoint);
        G4ThreeVector localDir = curDir;
        if (globalToLocal.IsRotated()) localDir = globalToLocal.TransformAxis(curDir);

        // Get the distance to leave this volume, in the given direction
        G4double arcLen = pCurVolume->GetLogicalVolume()->GetSolid()->DistanceToOut(localPos,localDir);

        //If we are already in a sensitive volume then simple calculation for distance
        if (pCurVolume->GetLogicalVolume()->GetSensitiveDetector())
        {
            //This calculates the distance along the direction of the track to the 
            //boundary of the current volume
            distToActArea = arcLen;
        }
        //Check for a pathological case
        else if (arcLen > maxStep)
        {
            distToActArea = -arcLen;
        }
        //Otherwise, we need to step until we find the next sensitive volume
        else
        {
            // Reset arc length to watch for embedded volumes
            arcLen = 0.;

            while(arcLen <= maxStep)
            {
                double safeStep     = 0.1; 
                double stepOverDist = 1000. * kCarTolerance;

                // Create the next point "just over the boundary"
                G4ThreeVector overPoint  = curPoint + (arcLen + stepOverDist) * curDir;

                //Use the G4 navigator to locate the current point, then compute the distance
                //to the current volume boundary
                pCurVolume = navigator->LocateGlobalPointAndSetup(overPoint, &curDir, true, true);

                // Make sure we have not stepped out of GLAST (how can this happen?)
                if (pCurVolume == 0) break;

                // If we have entered an active volume, then we are done
                if (pCurVolume->GetLogicalVolume()->GetSensitiveDetector())
                {
                    break;
                }

                // Get global to local transformation
                ////const G4AffineTransform& gblToLocal = navigator->GetGlobalToLocalTransform();
                ////localPos = gblToLocal.TransformPoint(overPoint);
                ////localDir = curDir;
                ////if (gblToLocal.IsRotated()) localDir = gblToLocal.TransformAxis(curDir);

                // Get the distance to leave this volume, in the given direction
                ////G4double distToOut = pCurVolume->GetLogicalVolume()->GetSolid()->DistanceToOut(localPos,localDir);

                // Otherwise, we need to step into the next volume and keep going
                double trackLen = navigator->ComputeStep(overPoint, curDir, maxStep, safeStep) + stepOverDist;

                //If we are right on the boundary then jump over and recalculate
                while(trackLen <= stepOverDist)
                {
                    // We are here because we are "on" the boundary (within the tolerance for that)
                    // Two reasons for the problem: 1) travelling almost exactly parallel to the boundary
                    // or, 2) not parallel
                    // Try to ascertain which here
                    G4bool              temp = false;
                    const G4ThreeVector localExitNrml = navigator->GetLocalExitNormal(&temp);
                    G4ThreeVector       trackDir = curDir;

                    if (globalToLocal.IsRotated()) trackDir = globalToLocal.TransformAxis(curDir);

                    double trkToExitAng = trackDir.dot(localExitNrml);

                    // Parallel track case
                    if (fabs(trkToExitAng) < kCarTolerance)
                    {
                        // One of the components of the direction vector is near zero (and below tolerance)
                        // and this is causing the step length calculation to screw up. The step length 
                        // calculation is done in local coordinates though...
                        //
                        // WHAT WE USED TO DO (THAT DIDN'T WORK ON LINUX)
                        // , so if we zero out the bad global
                        // component then the transformation will probably still result in a small but non-zero
                        // number... (and the G4 code tests on identically zero). 
                        // So... zero out the bad value in local coordinates
                        //if (fabs(trackDir.x()) < kCarTolerance) trackDir.setX(0.);
                        //if (fabs(trackDir.y()) < kCarTolerance) trackDir.setY(0.);
                        //if (fabs(trackDir.z()) < kCarTolerance) trackDir.setZ(0.);
            
                        // Now find the transformation from local back to global coordinates
                        //G4AffineTransform  localToGlobal = navigator->GetLocalToGlobalTransform();

                        // Update the current direction in global coordinates to this direction. 
                        // This will assure us that the value will be "correct" in the G4 step length calculation
                        //curDir = localToGlobal.TransformAxis(trackDir);

                        // Now what we do is to transform the localExitNrml vector to the global system
                        // and then shift the point by twice kCarTolerance away from the edge in the 
                        // opposite direction. This should leave us in the original volume, but now outside
                        // the tolerance.
                        // It's too bad the G4 code doesn't do this right...
                        G4AffineTransform localToGlobal = navigator->GetLocalToGlobalTransform();
                    
                        G4ThreeVector globalExitNrml = localToGlobal.TransformAxis(localExitNrml);

                        overPoint -= 2. * kCarTolerance * globalExitNrml;
                    }
                    // Surface tolerance case 
                    else 
                    {
                        stepOverDist += 1000. * kCarTolerance / fabs(trkToExitAng); // 10 * minimum tolerance in G4
                        overPoint     = curPoint + stepOverDist*curDir;
                    }

                    // Re-step with new points
                    pCurVolume = navigator->LocateGlobalPointAndSetup(overPoint, &curDir, true, true);
                    trackLen   = navigator->ComputeStep(overPoint, curDir, maxStep, safeStep) + stepOverDist;
                }

                // Test again, in the case that the above while loop was executed
                if (pCurVolume->GetLogicalVolume()->GetSensitiveDetector())
                {
                    break;
                }

                //Update distance stepped
                arcLen     += trackLen;
                curPoint    = curPoint + trackLen * curDir;
            }

            distToActArea = -arcLen;
        }
    }

    return distToActArea;
}

ParticleTransporter::ConstStepPtr ParticleTransporter::getStepAtArcLen(double arcLen) const
{
  // Purpose and Method: Goes through step list to find the step a distance arcLen from start
  // Inputs:  None
  // Outputs:  a constant iterator to the step
  // Dependencies: None
  // Restrictions and Caveats: None

    ConstStepPtr stepPtr = stepInfo.begin();
    double stepArcLen = 0.;

    // If arcLen is negative then want last step
    if (arcLen < 0.) stepPtr = stepInfo.end();
    // Otherwise, go through the steps until we pass the desired distance
    else
    {
        while(stepPtr < stepInfo.end())
        {
            double stepDist = (*stepPtr).GetArcLen();

            if (stepArcLen + stepDist > arcLen) break;

            stepArcLen += stepDist;
            stepPtr++;
        }
    }

    if(stepPtr!=stepInfo.begin()) --stepPtr;

    return stepPtr;
}


void ParticleTransporter::clearStepInfo()
{
  // Purpose and Method: Clears the vector of TransportStepInfo
  // objects. Typically before initiating a new step
  // Inputs:  None
  // Outputs:  None
  // Dependencies: None
  // Restrictions and Caveats: None

//  StepPtr stepPtr = stepInfo.begin();
//  while(stepPtr < stepInfo.end()) delete *stepPtr++;

  stepInfo.clear();

  return;
}


//Starting the the current volume, this builds the complete volume identifier
//string needed to specify where we are in the GLAST world
idents::VolumeIdentifier ParticleTransporter::constructId(const Hep3Vector& position, const Hep3Vector& direction, bool fudge) const
{
    // Purpose and Method: Constructs the complete volume indentifier from current
    //                     volume
    // Inputs: A pointer to a G4VPhysicalVolume object giving the current
    //         (starting) volume
    // Outputs:  a VolumeIdentifier
    // Dependencies: Requires that the volume - idents map has been found
    // Restrictions and Caveats: None

    G4Navigator*  navigator = m_TransportationManager->GetNavigatorForTracking();
    G4ThreeVector curPoint  = position;
    G4ThreeVector curDir    = direction;

    // If known to be on a boundary then back up to be sure to be inside the volume
    if (fudge)
    {
      G4ThreeVector displacement = -0.001 * curDir;
      curPoint += displacement; 
    }

    // Look up current volume and set tree above it
    G4VPhysicalVolume*  pCurVolume = navigator->LocateGlobalPointAndSetup(curPoint, 0, true, false);

    // Let's be sure that we are inside a valid volume (ie inside of GLAST - it can happen!)
    if (pCurVolume == 0) 
    {
        std::stringstream errorStream;
        errorStream << "ParticleTransporter::constructId given invalid initial conditions. pos: " 
                    << position << " dir: " << direction;
        throw std::domain_error(errorStream.str());
    }

    G4TouchableHistory* touchable  = navigator->CreateTouchableHistory();

    idents::VolumeIdentifier ret = constructId(touchable);

    delete touchable;

    return ret;
}

//Starting the the current volume, this builds the complete volume identifier
//string needed to specify where we are in the GLAST world
idents::VolumeIdentifier ParticleTransporter::constructId(G4TouchableHistory* theTouchable) const
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
    if (theTouchable)
	{
	    for( int i = 0; i<theTouchable->GetHistoryDepth() ; ++i) 
        {
	        const G4VPhysicalVolume* pVolume = theTouchable->GetVolume(i); 

            // Look up the identifier for this volume
            //VolumeIdentifier id = (*m_IdMap)[pVolume];
            VolumeIdentifier id = m_geometrySvc->getVolumeIdent(pVolume);

            // Add this volume's identifier to our total id
            ret.prepend(id);
	    }
	}

    return ret;
}
 
void ParticleTransporter::printStepInfo(std::ostream& str) const
{
    // Purpose and Method: Provides method for formatted outputting of tracking results
    // Inputs:  a reference to a standard output stream
    // Outputs:  None
    // Dependencies: None
    // Restrictions and Caveats: None


    str << "** Particle at: " << getStep(0).GetEntryPoint() << ", dir: " << getStep(0).GetDirection() << "\n";
    str << "   Step Count:  "  <<getNumberSteps() <<'\n';

    ConstStepPtr iter    = getStepStart();
    double       arcLen  = 0.;
    int          stepCnt = 0;
    double       x0sTot  = 0.;

    while(iter < getStepEnd())
    {
        TransportStepInfo        stepInfo   = *iter++;
        G4ThreeVector            position   = stepInfo.GetEntryPoint();
        G4ThreeVector            direction  = stepInfo.GetDirection();
        G4ThreeVector            entryPoint = position;
        const G4VPhysicalVolume* pCurVolume = stepInfo.GetVolume();
        G4Material*              pMaterial  = pCurVolume->GetLogicalVolume()->GetMaterial();
        G4String                 volName    = printVolName(pCurVolume);
        double                   matRadLen  = pMaterial->GetRadlen();
        double                   x0s        = 0.;
        double                   stepDist   = stepInfo.GetArcLen();

        if (matRadLen > 0.) x0s = stepDist / matRadLen;

        x0sTot += x0s;

        arcLen += stepDist;

        entryPoint -= stepDist * direction;

        //str << "   Step " << ++stepCnt << ", Volume: " << volName << "\n" 
        //    << ", position: " << position << ", step: " << stepDist << ", radlens: " << x0s << "\n";

        idents::VolumeIdentifier volId = constructId(position, direction);

        str << "   Step " << ++stepCnt << ", Volume: " << pCurVolume->GetName() << ":" << volId.name() 
             << ", entry position: " << entryPoint << "\n   step length: " << stepDist << ", radlens: " 
             << x0s << ", total: " << x0sTot << "\n";

        if (volId.size()!=9) continue;
        // check that it's really a TKR hit (probably overkill)
        if (!(volId[0]==0 && volId[3]==1 && volId[6] < 2)) continue; // !(LAT && TKR && Active)

        int trayNum = volId[4];
        int botTop  = volId[6];
        int view    = volId[5];
        int layer   = trayNum - 1 + botTop;

        str << "   ---> Sense Layer: " << layer << ", view: " << view << ", Bot/Top: " << botTop << "\n";          
    }
    str <<"   Total Step length: "<< arcLen <<"\n";

    return;
}

G4String ParticleTransporter::printVolName(const G4VPhysicalVolume* pCurVolume) const
{
    // Purpose and Method: Builds a string containing the entire volume name, from top 
    //                     mother volume to that specified by pCurVolume
    // Inputs:  Pointer to the lowest level G4VPhysicalVolume 
    // Outputs:  bool, returns true if target stop point reached, false otherwise
    // Dependencies: None
    // Restrictions and Caveats: None

    if (pCurVolume)
    {
    //    return "(" + pCurVolume->GetName() + printVolName(pCurVolume->GetMother()) + ")";
        return " ";
    }
    else return " ";
}

G4VPhysicalVolume* ParticleTransporter::getVolume(const Hep3Vector& position, const Hep3Vector& direction, bool fudge) const
{
    // Purpose and Method: Returns a pointer to the lowest level G4PhysicalVolume
    //                     which contains the point given in the variable position
    // Inputs:  The position and a bool to control use of a fudge to displace point
    //          off of the surface if at the edge of a volume.
    // Outputs:  bool, returns true if target stop point reached, false otherwise
    // Dependencies: Requires initialization via SetInitStep
    // Restrictions and Caveats: None

    G4Navigator*  navigator = m_TransportationManager->GetNavigatorForTracking();
    G4ThreeVector curPoint  = position;
    G4ThreeVector curDir    = direction;

    // If known to be on a boundary then back up to be sure to be inside the volume
    if (fudge)
    {
        G4ThreeVector displacement = -0.001 * curDir;
        curPoint += displacement; 
    }

    // Look up current volume and set tree above it
    G4VPhysicalVolume* pCurVolume = navigator->LocateGlobalPointAndSetup(curPoint, 0, true, false);

    return pCurVolume;
}
