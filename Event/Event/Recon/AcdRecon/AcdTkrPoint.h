#ifndef Event_ACDTKRPOINT_H
#define Event_ACDTKRPOINT_H

#include <vector>

#include "Event/TopLevel/Definitions.h"

#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

#include "idents/AcdId.h"
#include "geometry/Point.h"

#include "Event/Recon/TkrRecon/TkrTrackParams.h"

class MsgStream;

static const CLID& CLID_AcdTkrPointCol = InterfaceID("AcdTkrPointCol", 1, 0);

/**
*  @class Event::AcdTkrPoint
*
*  @brief TDS object which stores information where an extrapolated track exits the nominal ACD volume.  
*
*  This Point is calculated in 3D.  We include information the track paramterization at the point.
*  
*  The main access functions are:
*    - float getArcLength()
*      - which returns the arc-length along the track at which the point occurs.  
*        Postive values are given for intersections above the first track hit.
*    - const Point& point()
*      - which returns the intersection point in global coordinates
*    - int face() 
*      - which return the side of the ACD did we exited (0=Top, 1=-X, 2=-Y, 3=+X, 4=+Y, 5=Bottom)
*    - const Event::TkrTrackParams& paramsAtPoint()
*      - which returns the kalman propagated track parameters at the POINT
*  
*
*  \author Eric Charles
*
* $Header$
*/

namespace Event
{
  
  class AcdTkrPoint : virtual public ContainedObject {

  public:
    
    /// Default constructor.  
    AcdTkrPoint();
    
    /// Constructor for use in persistent -> transient conversion and reconstruction
    /// Takes arguements as they are stored in ROOT and caluclated by AcdReconAlg
    AcdTkrPoint(double arcLength, int face, 
		const Point& point, const Event::TkrTrackParams& paramsAtPoint);

    /// Destructor is trivial
    virtual ~AcdTkrPoint() {;}

    /// Return the arclength along the track at which the track cross the nominal ACD volume
    /// This is calculated in 3D.  
    inline float arcLength() const {
      return m_arcLength;
    }

    /// Return the point (in global coordinates)
    const Point& point() const {
      return m_point;
    }

    /// Which side of the ACD did we exited (0=Top, 1=-X, 2=-Y, 3=+X, 4=+Y, 5=Bottom)
    inline int face() const {
      return m_face;
    }

    /// Return the kalman propagated track parameters at the POINT
    const Event::TkrTrackParams& paramsAtPoint() const {
      return m_paramsAtPoint;
    }

    /// set all the values at once
    void set(double arcLength, int face, 
	     const Point& point, const Event::TkrTrackParams& paramsAtPoint);
    
    // set the individual values
    inline void setArcLength(float val) { m_arcLength = val; };
    inline void setFace(int val) { m_face = val; };
    inline void setPoint(const Point& val) { m_point = val; };
    inline void setParams(const Event::TkrTrackParams& val) { m_paramsAtPoint = val; };
    
    /// Print out this structure on a stream
    virtual void writeOut(MsgStream& stream) const;
    
  protected:

    /// Reset all the values to their defaults
    virtual void ini();
    
  private:
    
    /// The arclength at which this track crosses the nominal ACD volume
    float m_arcLength;

    /// which face of the ACD does the track leave by
    int m_face;
    
    /// The point of closest approach of the track to the tile or ribbon
    Point m_point;

    /// The track parameters (and covarience matrix) at the POINT
    Event::TkrTrackParams m_paramsAtPoint;
    
  };

   
  /*! 
   * @class AcdTkrPointCol
   *
   *  @brief TDS class a collection of AcdTkrPoint objects
   *  
   * It inherits from DataObject
   * (to be a part of Gaudi TDS) and from std::vector of pointers
   * to AcdTkrPoint objects. Some methods just rename the std::vector
   * methods with the same functionality for backward compartibility.
   *
   *
   * @author Eric Charles
   *
   * @todo replace this class by typedef to ObjectVector, make corresponding
   *       changes in AcdReconAlg
   */

  class AcdTkrPointCol : public DataObject, public std::vector<AcdTkrPoint*> 
  {
  public:
    
    /// Default constructor.  Builds empty collection
    AcdTkrPointCol() { clear();}
    
    /// "copy" constructor.  Take ownerships of a vector of AcdTkrPoint
    AcdTkrPointCol(const std::vector<AcdTkrPoint*>& acdhits);
  
    /// destructor - deleting the hits pointed
    /// by the vector elements
    ~AcdTkrPointCol() { delTkrPoints();}

    
    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_AcdTkrPointCol;}
    virtual const CLID& clID() const {return classID();}
        
    /// takes ownership of a vector AcdTkrPoint
    void init(std::vector<AcdTkrPoint*>& acdhits) {
      for ( std::vector<AcdTkrPoint*>::iterator itr = acdhits.begin(); itr != acdhits.end(); itr++ ) {
	push_back(*itr);
      }
    }

    /// add new AcdTkrPoint
    void add(AcdTkrPoint* cl) {push_back(cl);}
    
    /// get the number of points in collection
    int num()                  const {return size();}
    
    /// get pointer to the point with a given number 
    AcdTkrPoint * getTkrPoint(int i) const {return operator[](i);}
    
    /// delete all points pointed by the vector elements
    void delTkrPoints();
    
    /// write information for all points to the ascii file 
    /// for debugging purposes
    virtual void writeOut(MsgStream& stream) const;
    
  protected:
    
    /// does the same function as clear() 
    virtual void ini();
        
  };
  
}

#endif
