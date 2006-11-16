#ifndef GCRReconClasses_H
#define GCRReconClasses_H

#include "geometry/Point.h"
#include "geometry/Vector.h"

/**#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Vector/ThreeVector.h" 
#include "CLHEP/Vector/LorentzVector.h"
*/

#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "idents/CalXtalId.h"

#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/MsgStream.h"

//static const CLID& CLID_CalMipTrackVecCol = InterfaceID("CalMipTrackVecCol", 1, 0);

/**   
* @class GcrXtal
*
* 
*/

//-----------------------------------------------------------------------------------------------------------------
namespace Event 
{
    //-----------------------------------------------------------------------------------------------------------------
    // Define a GcrXtal class which will be used in GCRReconAlg
    class GcrXtal: virtual public ContainedObject 
    {
	private:
	    //Event::CalXtalRecData* m_xtalData;
	    idents::CalXtalId m_xtalId;
	    double                  m_pathLength;
	    

	    double m_closestFaceDist;
	    int m_crossedFaces;
	    Point m_entryPoint;
	    Point m_exitPoint;
	    

	public:
	    GcrXtal(){};

	    GcrXtal(idents::CalXtalId xtalId, double pathLength) : 
                	m_xtalId(xtalId), m_pathLength(pathLength), m_closestFaceDist(0.0), m_crossedFaces(0), m_entryPoint(Point()), m_exitPoint(Point()){};


	    GcrXtal(idents::CalXtalId xtalId, double pathLength,double closestFaceDist, int crossedFaces, Point entryPoint, Point exitPoint) : 
                	m_xtalId(xtalId), m_pathLength(pathLength), m_closestFaceDist(closestFaceDist), m_crossedFaces(crossedFaces), m_entryPoint(entryPoint), m_exitPoint(exitPoint) {};

	    ~GcrXtal() {};

	    void                   initialize(idents::CalXtalId xtalId, double pathLength, double closestFaceDist, int crossedFaces, Point entryPoint, Point exitPoint);

	    void                   setXtalId (idents::CalXtalId xtalId) {m_xtalId = xtalId ;}
	    void                   setPathLength  (double pathLength)         {m_pathLength  = pathLength;}
	    
	    void                   setClosestFaceDist (double closestFaceDist) {m_closestFaceDist = closestFaceDist ;}
	    void                   setCrossedFaces  (int crossedFaces)         {m_crossedFaces  = crossedFaces;}
	    void                   setEntryPoint (Point entryPoint) {m_entryPoint = entryPoint ;}
	    void                   setExitPoint  (Point exitPoint)         {m_exitPoint  = exitPoint;}



	    idents::CalXtalId getXtalId () const                         {return m_xtalId     ;}
	    double                  getPathLength () const               {return m_pathLength    ;}

	    double                  getClosestFaceDist () const          {return m_closestFaceDist    ;}
	    int                     getCrossedFaces () const             {return m_crossedFaces    ;}
	    Point                  getEntryPoint () const                {return m_entryPoint    ;}
	    Point                  getExitPoint () const                 {return m_exitPoint    ;}

	    /// Utilities 
	    void writeOut(MsgStream& log) const; 
	    std::ostream& fillStream( std::ostream& s ) const;
	    
	    void getReadableXedFaces(int xedFaces, int& inFace, int& outFace);
	    
    };
    

    // Define a vector of GcrXtals
    typedef std::vector<GcrXtal> GcrXtalVec;

    // Define a vector of GcrXtals
    typedef ObjectVector<GcrXtal>      GcrXtalCol;
    
    
    
    class GcrTrack:public DataObject
    {
    
    	private:
	
	    Vector m_direction;
	    Vector m_dirError;
	    Vector m_mcDirection;
	    Point m_calEntryPoint;
	    Point m_calExitPoint;



	public:
	    GcrTrack(){};

	    GcrTrack(Vector direction, Vector dirError, Point calEntryPoint, Point calExitPoint) : 
                	m_direction(direction), m_dirError(dirError), m_calEntryPoint(calEntryPoint), m_calExitPoint(calExitPoint) {};

	    GcrTrack(Vector direction, Vector mcDirection, Vector dirError, Point calEntryPoint, Point calExitPoint) : 
                	m_direction(direction), m_dirError(dirError), m_mcDirection(mcDirection), m_calEntryPoint(calEntryPoint), m_calExitPoint(calExitPoint) {};

	    ~GcrTrack() {};

	    void                   initialize(Vector direction, Vector dirError, Point calEntryPoint);

	    void                   setDirection (Vector direction)         {m_direction  = direction;}
	    void                   setDirError  (Vector dirError)         {m_dirError  = dirError;}
	    void                   setMcDir  (Vector mcDirection)         {m_mcDirection  = mcDirection;}
	    void                   setCalEntryPoint  (Point calEntryPoint)         {m_calEntryPoint  = calEntryPoint;}
	    void                   setCalExitPoint  (Point calExitPoint)         {m_calExitPoint  = calExitPoint;}
	   
	    Vector                 getDirection () const                   {return m_direction    ;}
	    Vector                 getDirError () const                    {return m_dirError    ;}
	    Vector                 getMcDir () const                       {return m_mcDirection    ;}
	    Point                  getCalEntryPoint () const               {return m_calEntryPoint    ;}
	    Point                  getCalExitPoint () const                {return m_calExitPoint    ;}

           
	    /// Utilities 
	    //void writeOut(MsgStream& log) const; 
	    //std::ostream& fillStream( std::ostream& s ) const;
    
    
    
    };
    
        // Define a vector of GcrXtals
    typedef ObjectVector<GcrTrack>      GcrTrackCol;

  
    
    
}
#endif

