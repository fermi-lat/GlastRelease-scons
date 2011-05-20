#ifndef Event_ACDTKRLOCALCOORDS_H
#define Event_ACDTKRLOCALCOORDS_H

#include "Event/TopLevel/Definitions.h"
#include "Event/Recon/TkrRecon/TkrTrackParams.h"
#include "idents/AcdId.h"

#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Geometry/Point3D.h"

class Point;

class MsgStream;

typedef HepGeom::Point3D<double> HepPoint3D;

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
    
    AcdTkrLocalCoords(int volume, float arcLength, 
                      const HepPoint3D& global);

    AcdTkrLocalCoords(int volume, float arcLength, float cosTheta, 
                      const HepPoint3D& global, 
                      const float localPosition[2], const float active[2], 
                      const CLHEP::HepSymMatrix& localCovProj, const CLHEP::HepSymMatrix& localCovProp);

    AcdTkrLocalCoords(float arcLength, float cosTheta, 
                      const HepPoint3D& global, 
                      const double localPosition[2], 
                      const CLHEP::HepSymMatrix& planeError);

    AcdTkrLocalCoords(const AcdTkrLocalCoords& other);

    virtual ~AcdTkrLocalCoords() {};
    
    AcdTkrLocalCoords& operator=(const AcdTkrLocalCoords& other);
      
    /// Direct access to parameters

    /// Which Volume of the tile or ribbon is this w.r.t.
    inline int   getLocalVolume()           const {return m_volume; }

    /// Location of hit in global coordinates
    inline const HepPoint3D& getGlobalPosition() const {return m_global; };

    /// Location (and errors) of hit in active distances
    inline float  getActiveX()               const {return m_active[0]; };
    inline float  getActiveY()               const {return m_active[1]; };

    /// Location (and errors) of hit in tile coordinates
    inline float  getLocalX()               const {return m_local[0]; };
    inline float  getLocalY()               const {return m_local[1]; };

    inline const CLHEP::HepSymMatrix& getLocalCovProj() const {return m_localCovProj; };
    inline const CLHEP::HepSymMatrix& getLocalCovProp() const {return m_localCovProp; };

    /// For backwards compatibility
    inline float getLocalXXCov() const { return m_localCovProj(1,1); }
    inline float getLocalXYCov() const { return m_localCovProj(1,2); }
    inline float getLocalYYCov() const { return m_localCovProj(2,2); }

    /// For backwards compatibility
    inline float getPathLength() const { return 0.; }

    /// Arclength along track to the detector plane
    inline float getArclengthToPlane()     const {return m_arcLengthPlane; }

    /// Angle of track relative to detector plane
    inline float  getCosTheta()            const {return m_cosTheta; }

    /// set everything at once
    void setLocalData(int volume, float arcLength, float cosTheta, 
                      const HepPoint3D& global, 
                      const float localPosition[2], const float active[2], 
                      const CLHEP::HepSymMatrix& localCovProj, const CLHEP::HepSymMatrix& localCovProp);
    
    /// set everything at once, old version
    void setLocalData(const float localPosition[2],
                      float pathLength, float cosTheta, 
                      int region, const CLHEP::HepSymMatrix& planeError);
    
    /// set stuff from old version of AcdTkrPoint
    void setLocalData(float arcLength, int face, const Point& point, const Event::TkrTrackParams& params);

    /// copy in everything at once
    void copy(const AcdTkrLocalCoords& other);
    
    virtual void writeOut(MsgStream& stream) const;
    
  protected:

    virtual void ini();
    
  private:
    
    ///  Which volume of the Tile or Ribbon does this occur in 
    int      m_volume;
    
    ///  Arclength to plane
    float    m_arcLengthPlane;

    ///  Angle of track w.r.t. detector plane
    float    m_cosTheta;

    /// Global position of expected hit 
    HepPoint3D m_global;

    ///  Position of the hit in Local Coords
    float    m_local[2];

    ///  Active distance of the hit 
    float    m_active[2];

    ///  Covariance terms in expected LocalCoords
    CLHEP::HepSymMatrix m_localCovProj;
 
        ///  Covariance terms in expected intersection
    CLHEP::HepSymMatrix m_localCovProp;

  };

}

#endif
