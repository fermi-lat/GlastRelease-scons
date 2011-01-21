#ifndef EVENT_ACDPOCADATA_H
#define EVENT_ACDPOCADATA_H

/** 
 * @class Event::AcdPocaData
 *
 * @brief TDS object information about the Point of Closest Approach (POCA) between an extrapolated track
 *  and the Acd (either a detector element or a gap).  
 *
 *  This POCA is calculated in 3D.  The doca is defined to be positive if the track goes inside the active distance 
 *  or gap and negative otherwise.  Both AcdTkrHitPoca and AcdTkrGapPoca inherit from this class.
 *
 *  The main access functions are:
 *    - float getArcLength()
 *      - which returns the arc-length along the track at which the POCA occurs.  
 *        Postive values are given for intersections above the first track hit.
 *    - float getDoca(), float getDocaErr()
 *      - which returns the distance of closest approach.  This is usually expressed as an active distance
 *        Positive values mean the track went into the volume in question, 
 *        Negative values mean it missed the volume.
 *        This magnitude is the distance to the edge of the volume.
 *    - const Point& getPoca()
 *      - which returns the Point of Closest Approach (POCA) in global coordinates
 *    - const Vector& getPocaVector()
 *      - which returns the  vector between the POCA and the closest edge of the volume in question
 *
 * @author Eric Charles
 *
 * $Header$
 */

#include "geometry/Vector.h"
#include "geometry/Point.h"

class MsgStream;

namespace Event {

  class AcdPocaData {
    
  public:
    
    /// Default constructor.  
    AcdPocaData();
    
    /// Constructor for use in transient -> persistent conversion 
    /// Takes arguements as they are stored in ROOT
    AcdPocaData(int volume, int region, float arcLength, 
                float doca, float docaErrProj, float docaErrProp,
                const Point& poca, const Vector& voca);
    
    /// Copy constructor
    AcdPocaData(const AcdPocaData& other);
    
    virtual ~AcdPocaData() {;}
    
    /// Assignment operator
    virtual AcdPocaData& operator=(const AcdPocaData& other);
    
    /// Return the sub-volume of the tile or ribbon associated with this poca
    inline int getVolume() const {
      return m_volume;
    }

    /// Region the region of the sub-volume associated with the poca
    inline int getRegion() const {
      return m_region;
    }

    /// Return the arcLength along the track at which the POCA occurs
    /// This is calculated in 3D.  See AcdPocaTool for details  
    inline float getArcLength() const {
      return m_arcLength;
    }
    
    /// Return the distance of closest approach (doca)
    /// This is calculated in 3D.  See AcdPocaTool for details  
    inline float getDoca() const {
      return m_doca;
    }
    
    /// Return the distance of closest approach (doca)
    /// This is calculated in 3D.  See AcdPocaTool for details  
    inline float getDocaErr() const {
      return m_docaErr_proj;
    }

    /// Return the error on distance of closest approach (doca)
    inline float getDocaErrProj() const {
      return m_docaErr_proj;
    }

    /// Return the error on distance of closest approach (doca)
    /// This is calulated the full Kalman Filter propagation of the track
    /// to the POCA and the projection of the propagated covarience matrix 
    /// along the line of the doca.
    inline float getDocaErrProp() const {
      return m_docaErr_prop;
    }
    
    /// Return the POCA (in global coordinates)
    inline const Point& getPoca() const {
      return m_poca;
    }
    
    /// Return vector between the POCA and the closest edge of the detector element or GAP 
    inline const Vector& getPocaVector() const {
      return m_voca;
    }
    
    
    /// set all the values
    void setPocaData(int volume, int region, float arcLength, 
             float doca, float docaErrProj, float docaErrProp,
             const Point& poca, const Vector& pocaVector);

    /// set all the values, old version
    void setPocaData(float arcLength, float doca, float docaErr, 
             const Point& poca, const Vector& pocaVector);
    
    void setPocaData(const AcdPocaData& other);
    
    /// set individaul values
    
    virtual void writeOut(MsgStream& stream ) const;
    
  protected:

    /// reset all the values to their default (lowercase to avoid conflict with TObject::Clear())
    void ini();
    
  private:
    
    /// Which volume of the Tile or Ribbon does this occur w.r..t
    int   m_volume;

    /// Which region of the poca occur w.r.t. (ie, which edge or corner of the tile or ribbon)
    int   m_region;

    /// The arclength at which the poca occurs
    float m_arcLength;
    
    /// The Distance of Closest Approach between track and tile or ribbon
    float m_doca;                 
    
    /// The Error on the DOCA, this is the track error ellipsoid project along the doca direction
    float m_docaErr_proj;

    /// The Error on the DOCA, this is the track error ellipsoid project along the doca direction
    float m_docaErr_prop;
    
    /// The Point of Closest Approach
    Point m_poca;
    
    /// The vector from the POCA to the Detector Element edge
    Vector m_voca;
    
  };
    
};
#endif
