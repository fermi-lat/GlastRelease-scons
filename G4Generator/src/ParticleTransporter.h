#ifndef ParticleTransporter_h
#define ParticleTransporter_h 1

#ifdef WIN32 // please leave here; THB 
#include <float.h>
#endif

#include "G4SteppingManager.hh"
#include "DetectorConstruction.h"
#include "TransportStepInfo.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"

#include "globals.hh"
#include "g4std/vector"

#include "GlastSvc/Reco/IKalmanParticle.h"

#include <vector>

/** 
* @class ParticleTransporter
*
* @brief Geant4 class for particle transport
*
* This class is the implementation of an algorithm for tracking a particle from
* a starting point to a finish point defined by either a desired arc length, or
* by reaching (and stepping through) the "next" sensitive volume. It keeps track
* of the results of each step in a vector of TransportStepInfo objects.
*
* @author Tracy Usher
*
*/

class ParticleTransporter  
{
public: 
    ParticleTransporter();
   ~ParticleTransporter();

   /**  
    * Sets the starting point and direction and will do an initial step of
    * distance step @param start - the initial starting point @param dir - the
    * initial direction @param step - the initial step length to take
    */
    void  setInitStep(const Point& start, const Vector& dir, 
                      const double step);

    /// Performs the actual stepping
    bool               transport();

    /// Methods for retrieving information after tracking
    int                getNumberSteps()  const {return stepInfo.size();}
    ConstStepPtr       getStepStart()    const {return stepInfo.begin();}
    ConstStepPtr       getStepEnd()      const {return stepInfo.end();}

    TransportStepInfo* getLastStep()     const {return stepInfo.back();}
    //for now
    TransportStepInfo* getPrevBoundary() const {return stepInfo.back();} 

    HepPoint3D         getStartPoint()   const {return startPoint;}
    HepVector3D        getStartDir()     const {return startDir;}
    double             getMaxArcLen()    const {return maxArcLen;}

    void printStepInfo(std::ostream& str=std::cout ) const;  
private:
    //Private methods
    //Determines the distance along the track to the edge of start volume
    double     minStepSize(G4SteppingManager* stepManager);
    //Make sure list is cleared when tracking started
    void       clearStepInfo();
    //Use this to in printStepInfo to print volume list name
    G4String   printVolName(const G4VPhysicalVolume* pCurVolume) const;

    //Private Data
    //Vector keeping track of the individual step information
    StepVector stepInfo;

    //Initial start point, direction and desired step size
    HepPoint3D  startPoint;
    HepVector3D startDir;
    double      maxArcLen;
};

#endif




