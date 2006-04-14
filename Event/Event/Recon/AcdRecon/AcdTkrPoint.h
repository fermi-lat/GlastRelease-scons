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
*  @class AcdTkrPoint
*
*
*  @brief This class stores information about the Point of Closest Approach (POINT) between an extrapolated track
*  and a hit Acd element (tile or ribbon).  This POINT is calculated in 3D.  The doca is defined to be positive 
*  if the track goes inside the active distance and negative otherwise
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

    /// Return the POINT (in global coordinates)
    const Point& point() const {
      return m_point;
    }

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
   *  @brief TDS class  to store the results of the reconstruction performed
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
   *       changes in AcdTrkIntersectionAlg
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
