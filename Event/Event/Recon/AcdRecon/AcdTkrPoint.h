#ifndef Event_ACDTKRPOINT_H
#define Event_ACDTKRPOINT_H

#include "Event/Recon/AcdRecon/AcdTkrLocalData.h"
#include "Event/Recon/TkrRecon/TkrTrackParams.h"

#include <vector>

#include "Event/TopLevel/Definitions.h"

#include "Event/Recon/TkrRecon/TkrTrackParams.h"

#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

#include "idents/AcdId.h"
#include "geometry/Point.h"


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
  
  class AcdTkrPoint : virtual public ContainedObject, public AcdTkrLocalCoords {

  public:
    
    /// Default constructor.  
    AcdTkrPoint();
    
    /// Copy c'tor
    AcdTkrPoint(const AcdTkrPoint& other);

    /// Constructor for use in persistent -> transient conversion and reconstruction
    /// Takes arguements as they are stored in ROOT and caluclated by AcdReconAlg
    AcdTkrPoint( int trackIndex,
                 int volume, float arcLength, float cosTheta, 
                 const HepPoint3D& global, const float localPosition[2], 
                 const CLHEP::HepSymMatrix& localCovProj, const CLHEP::HepSymMatrix& localCovProp);

    /// Old Constructor for backwards compatiblity
    AcdTkrPoint( float arcLength, int volume,
                 const HepPoint3D& global, const Event::TkrTrackParams& params );
    
    /// Destructor is trivial
    virtual ~AcdTkrPoint() {;}

    /// Assignment operator
    Event::AcdTkrPoint& operator=(const Event::AcdTkrPoint& other);

    /// Return the track index
    inline int getTrackIndex() const {
      return m_trackIndex;
    }

    /// for backwards compatibility
    Point point() const { 
      return Point(AcdTkrLocalCoords::getGlobalPosition().x(), 
                   AcdTkrLocalCoords::getGlobalPosition().y(),
                   AcdTkrLocalCoords::getGlobalPosition().z()); }

    /// for backwards compatibility
    const Event::TkrTrackParams& paramsAtPoint() const {
      static const Event::TkrTrackParams nullParams;
      return nullParams;
    }

    /// for backwards compatibility
    inline float arcLength() const {
      return AcdTkrLocalCoords::getArclengthToPlane();
    }

    /// for backwards compatibility
    inline int face() const {
      return AcdTkrLocalCoords::getLocalVolume();
    }
 

    /// set all the values at once
    void set( int trackIndex );
    
    /// Print out this structure on a stream
    virtual void writeOut(MsgStream& stream) const;
    
  protected:

    /// Reset all the values to their defaults
    virtual void ini();
    
  private:
    
    /// Which Track
    int   m_trackIndex;
    
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
