#ifndef ParticleTransporter_h
#define ParticleTransporter_h 1

#ifdef WIN32 // please leave here; THB 
#include <float.h>
#endif

//Necessary G4 stuff
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"

//Get the geometry information
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "idents/VolumeIdentifier.h"
#include "G4Generator/IG4GeometrySvc.h"

#include <vector>

/** 
* @class ParticleTransporter
*
* @brief Geant4 class for particle transport
*
* This class implements an algorithm for tracking a particle from a starting point
* to a finish point, using the Geant4 version of the GLAST detmodel geometry and 
* the Geant4 volume routines for doing the tracking. It keeps track of the exit 
* positions and arclength traversed through each volume it encounters, and then 
* provides methods for recovering specific pieces of information for those volumes,
* or for the total track. 
*
* @author Tracy Usher
*
*/

class ParticleTransporter  
{
public: 
    ParticleTransporter(const G4TransportationManager* TransportationManager, IG4GeometrySvc* geoSvc);
   ~ParticleTransporter();

    /// @brief Internal class for storing information for each step taken
    class TransportStepInfo
    {
    public:
        TransportStepInfo(const G4ThreeVector& coords, G4double step) :
                          position(coords),arcLen(step) {};
       ~TransportStepInfo() {};

        G4ThreeVector       GetCoords() const {return position;}
        G4double            GetArcLen() const {return arcLen;}

        void SetCoords(const G4ThreeVector& coords) {position   = coords;}
        void SetArcLen(const G4double newArcLen)    {arcLen     = newArcLen;}
    
    private:
        G4ThreeVector      position;
        G4double           arcLen;
    };

    /// Useful typedefs for using the above internal class
    typedef std::vector<TransportStepInfo>                 StepVector;
    typedef std::vector<TransportStepInfo>::iterator       StepPtr;
    typedef std::vector<TransportStepInfo>::const_iterator ConstStepPtr;

    /// @brief Initializes the starting point (@param start) and direction (@param dir)
    /// before attempting to transport a track through the geometry.
    void  setInitStep(const Point& start, const Vector& dir);

    /// @brief Transports (steps) the track through an arclength @param step
    bool                     transport(const double step = -1.);

    /// @brief Methods for returning information on the stepping through the volumes
    int                      getNumberSteps()               const {return stepInfo.size();}
    /// Iterators
    ConstStepPtr             getStepStart()                 const {return stepInfo.begin();}
    ConstStepPtr             getStepEnd()                   const {return stepInfo.end();}
    ConstStepPtr             getStepAtArcLen(double arcLen) const;
    /// Objects
    TransportStepInfo        getStep(int stepIdx)           const {return stepInfo[stepIdx];}
    TransportStepInfo        getLastStep()                  const {return stepInfo.back();}
    TransportStepInfo        getPrevBoundary()              const {return stepInfo.back();} 

    /// @brief Returns the current volume (include tree above it) at given position. 
    /// If position is on a volume boundary, @param fudge = true will cause it to be
    /// moved slightly off the boundary.
    G4VPhysicalVolume*       getVolume(const Hep3Vector& position, bool fudge=false) const;

    /// @brief Given a volume, construct the associated volume id
    idents::VolumeIdentifier constructId(const Hep3Vector& position, bool fudge=false) const;

    /// @brief Return starting point information
    Point                    getStartPoint()                const {return startPoint;}
    Vector                   getStartDir()                  const {return startDir;}

    /// @brief Methods to return information after the last step taken
    double                   getTotalArcLen()               const;
    double                   insideActiveArea()             const;
    double                   insideActiveLocalX()           const;
    double                   insideActiveLocalY()           const;

    /// @brief Provide ability to print out information from each volume
    void                     printStepInfo(std::ostream& str=std::cout ) const;    

    /// @brief Methods to intialize the class
    void                     setTransportationManager(const G4TransportationManager* TransportationManager)
                             {m_TransportationManager = TransportationManager;}
    void                     setGeometrySvc(IG4GeometrySvc* geoSvc) {m_geometrySvc = geoSvc;}
private:
    //Private methods
    //Determines the distance along the track to the edge of start volume
    double                   minStepSize(const Point& start, const Vector& dir) const;
    //Make sure list is cleared when tracking started
    void                     clearStepInfo();
    //This function will step a given arc length
    bool                     StepAnArcLength(const double arcLen);
    //This function will step to the next "sensitive" plane
    bool                     StepToNextPlane();
    //Use this to in printStepInfo to print volume list name
    G4String                 printVolName(const G4VPhysicalVolume* pCurVolume) const;
    //This is used to find an "SiLadders" volume in the nested heirarchy
    const G4VPhysicalVolume* findSiLadders(G4TouchableHistory* theTouchable) const;
    //This used to find the distance to nearest edge in a given direction
    double                   distanceToEdge(const Vector& dir) const;
    //Gives the distance to the closest edge in a given and opposite direction
    double                   distanceToClosestEdge(const Vector& dir) const;
    /// brief Given a volume, construct the associated volume id
    idents::VolumeIdentifier constructId(G4TouchableHistory* theTouchable) const;


    //Private Data
    //Vector keeping track of the individual step information
    StepVector stepInfo;

    //Initial start point, direction and desired step size
    Point      startPoint;
    Vector     startDir;

    /// This is needed for associating Geant4 volumes to Glast identifiers
    IG4GeometrySvc* m_geometrySvc;

	/// Here we maintain a pointer to the G4 Transportation Manager
	const G4TransportationManager* m_TransportationManager;
};

#endif




