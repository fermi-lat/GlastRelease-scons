#ifndef Event_ACDTKRINTERSECTION_H
#define Event_ACDTKRINTERSECTION_H

#include <vector>
#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "Event/TopLevel/Definitions.h"
#include "idents/AcdId.h"

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

class MsgStream;
namespace CLHEP {
class HepMatrix;
};

static const CLID& CLID_AcdTkrIntersectionCol = InterfaceID("AcdTkrIntersectionCol", 1, 0);

/**
*  @class Event::AcdTkrIntersection
*
*  @brief The class stores information about where a track extrapolation crosses an element
*  of the GEANT detector mode.  
*
*  This class is mainly useful for measuring ACD instrument performance since it doesn't 
*  depend on the ACD so much as on the tracking.  Ie. AcdTkrIntersections are made even
*  if the ACD doesn't fire.
*
*  The main access functions are:
*    - const idents::AcdId& getTileId()  
*      - which returns the ID of the hit element
*    - const int getTrackIndex()  
*      - which returns the index of the track which did the hitting
*    - Point& getGlobalPosition() 
*      - which return the location of the hit in global coordinates
*    - double getLocalX(),double getLocalY()  
*      - which return the location of the hit in local cooridanates
*    - double getLocalXXCov() , double getLocalXYCov() , double getLocalYYCov()
*      - which return the terms of the covarience matrix projected into the place of the ACD element
*    - double getArcLengthToIntersection()
*      - which returns the arc-length along the track at which the hit occurs. 
*        Postive values are given for intersections above the first track hit
*    - double getPathLengthInTile()
*      - which returns the pathlength of track in detector element
*    - double getCosTheta()
*      - which returns the angle of the track relative to detector plane
*    - unsigned char tileHit()
*      - which returns a mask to say if the tile was hit
*        -AcceptMap_pmtA = 0x01
*        -AcceptMap_pmtB = 0x01
*        -Veto_pmtA      = 0x04
*        -Veto_pmtB      = 0x08
*        -CNO_pmtA       = 0x10
*        -CNO_pmtB       = 0x20
*  
*  \author Eric Charles
*
* $Header$
*/

namespace Event
{
  
  class AcdTkrIntersection 
  {
    
  public:

    AcdTkrIntersection();
    
    AcdTkrIntersection(const idents::AcdId& acdId, int trackIndex, 
                       const Point& globalPosition, 
               const double localPosition[2], const CLHEP::HepMatrix& localCovMatrix,
                       double arcLengthToIntersection, double pathLengthInTile,
                       unsigned char tileHit, double cosTheta);

    virtual ~AcdTkrIntersection() {};
    
    /// Direct access to parameters

    /// Which tile should be hit
    inline const idents::AcdId& getTileId()      const {return m_tileId; };
    /// Which track did the hitting
    inline int     getTrackIndex()          const {return m_trackIndex; };

    /// Location of hit in global coordinates
    inline const Point& getGlobalPosition()      const {return m_location;    };
    inline Point&       getGlobalPosition()            {return m_location;    };
    
    /// Location (and errors) of hit in tile coordinates
    inline double  getLocalX()              const {return m_localX; };
    inline double  getLocalY()              const {return m_localY; };
    inline double  getLocalXXCov()          const {return m_localXXCov; };
    inline double  getLocalYYCov()          const {return m_localYYCov; };
    inline double  getLocalXYCov()          const {return m_localXYCov; };

    /// Distance along track from first hit to tile intersection
    inline double  getArcLengthToIntersection() const { return m_arcLengthToIntersection; } ;
    /// Path length of track through tile
    inline double  getPathLengthInTile()    const { return m_pathlengthInTile; } ;
    /// Angle of track w.r.t. detector element
    inline double getCosTheta() const { return m_cosTheta; }

    /// mask to say if the tile was hit
    inline unsigned char tileHit() const { return m_tileHit; };
    
    /// set everything at once
    void set(const idents::AcdId& acdId, int trackIndex, 
             const Point& globalPosition, 
             const double localPosition[2], const HepMatrix& localCovMatrix,
             double arcLengthToIntersection, double pathLengthInTile,
             unsigned char tileHit, double cosTheta);

    virtual void writeOut(MsgStream& stream) const;
    
  protected:

    virtual void ini();
    
  private:
    
    /// ID of hit tile
    idents::AcdId m_tileId;

    /// index of the related track
    int       m_trackIndex;

    /// Global position of expected hit 
    Point     m_location;
    ///  X Position of the expected hit in Tile Coords
    double    m_localX;
    ///  Y Position of the expected hit in Tile Coords
    double    m_localY;
    
    ///  Covariance terms in expected intersection
    double    m_localXXCov;        // local X Error squared  (x for Top, x for +-y planes, y for +- x plane) 
    double    m_localYYCov;        // local Y Error squared  (y for Top, z for +-x and +-y planes)
    double    m_localXYCov;        // correlation term of local X-Y error   

    ///  Distance from first hit to intersection
    double    m_arcLengthToIntersection;
  
    ///  Pathlength through the ACD tile
    double    m_pathlengthInTile;
    
    ///  Mask to store tile hit
    unsigned char m_tileHit;

    /// Angle of track w.r.t. detector element
    double    m_cosTheta;

  };

   
  /*! 
   * @class AcdTkrIntersectionCol
   *
   *  @brief TDS class a collection of AcdTkrIntersection objects
   *  
   * It inherits from DataObject
   * (to be a part of Gaudi TDS) and from std::vector of pointers
   * to AcdTkrIntersection objects. Some methods just rename the std::vector
   * methods with the same functionality for backward compartibility.
   *
   *
   * @author Eric Charles
   *
   * @todo replace this class by typedef to ObjectVector, make corresponding
   *       changes in AcdReconAlg
   */
    
    
  class AcdTkrIntersectionCol : public DataObject, public std::vector<AcdTkrIntersection*> 
  {
  public:

    /// Fill from a vector of intersections, essentially a copy c'tor
    AcdTkrIntersectionCol(const std::vector<AcdTkrIntersection*>& acdTkrIntersections);

    /// Null c'tor
    AcdTkrIntersectionCol() { clear();}
    
    /// destructor - deleting the clusters pointed
    /// by the vector elements
    ~AcdTkrIntersectionCol() { delIntersections();}
        
    
    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_AcdTkrIntersectionCol;}
    virtual const CLID& clID() const {return classID();}
    
    /// add new cluster
    void add(AcdTkrIntersection* cl) {push_back(cl);}
    
    /// get the number of clusters in collection
    int num()                  const {return size();}
    
    /// get pointer to the cluster with given number 
    AcdTkrIntersection* getIntersection(int i) const {return operator[](i);}
    
    /// delete all intersections pointed by the vector elements
    void delIntersections();
    
    /// write information for all clusters to the ascii file 
    /// for debugging purposes
    virtual void writeOut(MsgStream& stream) const;
    
  protected:
    
    /// does the same function as clear() 
    virtual void ini();
        
  };
    
};

#endif
