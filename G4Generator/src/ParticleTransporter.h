#ifndef ParticleTransporter_h
#define ParticleTransporter_h 1

#ifdef WIN32 // please leave here; THB 
#include <float.h>
#endif

//Class def for keeping track of step info
#include "TransportStepInfo.h"

//Necessary G4 stuff
#include "G4Navigator.hh"
#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include "g4std/vector"

//Get the geometry information
#include "geometry/Point.h"
#include "geometry/Vector.h"

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
    int                getNumberSteps()   const {return stepInfo.size();}
    ConstStepPtr       getStepStart()     const {return stepInfo.begin();}
    ConstStepPtr       getStepEnd()       const {return stepInfo.end();}

    TransportStepInfo* getLastStep()      const {return stepInfo.back();}
    //for now
    TransportStepInfo* getPrevBoundary()  const {return stepInfo.back();} 

    Point              getStartPoint()    const {return startPoint;}
    Vector             getStartDir()      const {return startDir;}
    double             getMaxArcLen()     const {return maxArcLen;}
    double             insideActiveArea() const;
    double             insideActiveLocalX() const;
    double             insideActiveLocalY() const;

    void               printStepInfo(std::ostream& str=std::cout ) const;    
private:
    //Private methods
    //Determines the distance along the track to the edge of start volume
    double             minStepSize(const Point& start, const Vector& dir) const;
    //Make sure list is cleared when tracking started
    void               clearStepInfo();
    //This function will step a given arc length
    bool               StepAnArcLength();
    //This function will step to the next "sensitive" plane
    bool               StepToNextPlane();
    //Use this to in printStepInfo to print volume list name
    G4String           printVolName(const G4VPhysicalVolume* pCurVolume) const;
    //This is used to find an "SiLadders" volume in the nested heirarchy
    G4VPhysicalVolume* findSiLadders(const G4TouchableHistory* pHistory) const;
    //This used to find the distance to nearest edge in a given direction
    double             distanceToEdge(const Vector& dir) const;
    //Gives the distance to the closest edge in a given and opposite direction
    double             distanceToClosestEdge(const Vector& dir) const;


    //Private Data
    //Vector keeping track of the individual step information
    StepVector stepInfo;

    //Initial start point, direction and desired step size
    Point      startPoint;
    Vector     startDir;
    double     maxArcLen;
    double     minimumStep;
};

#endif




