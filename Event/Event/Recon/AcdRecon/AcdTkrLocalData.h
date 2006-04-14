#ifndef Event_ACDTKRLOCALCOORDS_H
#define Event_ACDTKRLOCALCOORDS_H

#include <vector>
#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "Event/TopLevel/Definitions.h"
#include "idents/AcdId.h"

#include "CLHEP/Matrix/Matrix.h"

class MsgStream;


/**
*  @class AcdTkrLocalCoords
*
*
*  @brief This class stores the results of the reconstruction performed
*  in AcdTkrIntersectAlg. It contains the reconstructed data for one
*  expected LocalCoords of a track with the calorimeter. 
*  
*  \author Eric Charles
*
* $Header$
*/

namespace Event
{
  
  class AcdTkrLocalCoords 
  {

  public:
    
    AcdTkrLocalCoords();
    
    AcdTkrLocalCoords(const float localPosition[2], float pathLength, float cosTheta, 
		      int region, const CLHEP::HepMatrix& localCovMatrix);

    AcdTkrLocalCoords(const AcdTkrLocalCoords& other);

    virtual ~AcdTkrLocalCoords() {};
    
    AcdTkrLocalCoords& operator=(const AcdTkrLocalCoords& other);
      
    /// Direct access to parameters

    /// Location (and errors) of hit in tile coordinates
    inline float  getActiveX()             const {return m_activeX; };
    inline float  getActiveY()             const {return m_activeY; };
    inline float  getLocalXXCov()          const {return m_localXXCov; };
    inline float  getLocalYYCov()          const {return m_localYYCov; };
    inline float  getLocalXYCov()          const {return m_localXYCov; };

    /// Angle of track relative to detector plane
    inline float  getCosTheta()            const {return m_cosTheta; }

    /// Pathlength of track in detector element
    inline float  getPathLength()          const {return m_pathLength; }

    /// A code which tells which region of the tile was hit
    inline int    getRegion()              const {return m_region; }

    /// set everything at once
    void set(const float localPosition[2], float pathLength, float cosTheta, 
	     int region, const CLHEP::HepMatrix& localCovMatrix);

    /// copy in everything at once
    void copy(const AcdTkrLocalCoords& other);
    
    // set the individual values (uncomment as needed)
    //inline void setLocalPosition(const float localPosition[2]) {
    //  m_localX = localPosition[0];
    //  m_localY = localPosition[1];
    //};
    //inline void setRegion(int val) {m_region = val;}
    //inline void setCosTheta(float val) { m_cosTheta = val; }
    //inline void setPathLength(float val) { m_pathLength = val; }
    inline void setCovTerms(float XX, float YY, float XY) {
      m_localXXCov = XX;
      m_localYYCov = YY;
      m_localXYCov = XY;      
    }

    virtual void writeOut(MsgStream& stream) const;
    
  protected:

    virtual void ini();
    
  private:
    
    ///  X Position of the expected hit in Local Coords
    float    m_activeX;
    ///  Y Position of the expected hit in Local Coords
    float    m_activeY;
    
    ///  Pathlength of track through tile
    float    m_pathLength;

    ///  Angle of track w.r.t. detector plane
    float    m_cosTheta;

    ///  Code that tells which part of tile was hit
    int      m_region;    

    ///  Covariance terms in expected LocalCoords
    float    m_localXXCov;        // local X Error squared  (x for Top, x for +-y planes, y for +- x plane) 
    float    m_localYYCov;        // local Y Error squared  (y for Top, z for +-x and +-y planes)
    float    m_localXYCov;        // correlation term of local X-Y error   

  };

}

#endif
