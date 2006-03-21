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
*  @class AcdTkrIntersection
*
*
*  @brief This class stores the results of the reconstruction performed
*  in AcdTkrIntersectAlg. It contains the reconstructed data for one
*  expected intersection of a track with the calorimeter. 
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
    
    AcdTkrIntersection(const idents::AcdId& acdId, int trackIndex, 
		       const Point& globalPosition, 
               const double localPosition[2], const CLHEP::HepMatrix& localCovMatrix,
		       double arcLengthToIntersection, double pathLengthInTile,
		       unsigned char tileHit);

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

    /// mask to say if the tile was hit
    inline unsigned char tileHit() const { return m_tileHit; };
    
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

  };

   
  /*! 
   * @class AcdTkrIntersectionCol
   *
   *  @brief TDS class  to store the results of the reconstruction performed
   *  in AcdTkrIntersectAlg
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
   *       changes in AcdTrkIntersectionAlg
   */
    

  class AcdTkrIntersectionCol : public DataObject, public std::vector<AcdTkrIntersection*> 
  {
  public:
    
    AcdTkrIntersectionCol(const std::vector<AcdTkrIntersection*>& acdTkrIntersections);

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

}

#endif
