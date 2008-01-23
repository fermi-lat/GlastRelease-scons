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
*  @class Event::AcdTkrLocalCoords
*
*  @brief TDS object which stores information about where a track crossed the plane of an ACD element.
*
*  This information is given in 2D.  The active distances are reported in the ACD element frame.
*  Both AcdTkrHitPoca and AcdTkrGapPoca inherit from this class.
*
*  The main access functions are:
*    - float getActiveX() , float  getActiveY()
*      - which return the active distance int the local frame of the ACD element
*        Positive values mean the track went into the volume in question, 
*        Negative values mean it missed the volume.
*        This magnitude is the distance to the edge of the volume.
*    - float getLocalXXCov() , float getLocalXYCov() , float getLocalYYCov()
*      - which return the terms of the covarience matrix projected into the place of the ACD element.
*    - float getCosTheta()
*      - which returns the angle of the track relative to detector plane
*    - float getPathLength()
*      - which returns the pathlength of track in detector element
*    - int getRegion()
*      - which returns a code telling which region of the tile was hit
*
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

    /// set the covarience terms
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
