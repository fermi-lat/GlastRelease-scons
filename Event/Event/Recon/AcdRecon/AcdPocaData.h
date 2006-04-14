#ifndef EVENT_ACDPOCADATA_H
#define EVENT_ACDPOCADATA_H

/** 
 * @class AcdPocaData
 * @brief TDS object information about the Point of Closest Approach (POCA) between an extrapolated track
 *  and the Acd (either a detector element or a gap).  This POCA is calculated in 3D.  
 *  The doca is defined to be positive if the track goes inside the active distance or gap
 *  and negative otherwise
 *
 * This class should be a duplicate of Event::AcdPocaData
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
    AcdPocaData(float arcLength, float doca, float docaErr, 
		const Point& poca, const Vector& pocaVector);
    
    /// Copy constructor
    AcdPocaData(const AcdPocaData& other);
    
    virtual ~AcdPocaData() {;}
    
    /// Assignment operator
    virtual AcdPocaData& operator=(const AcdPocaData& other);
    
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
    
    /// Return the error on distance of closest approach (doca)
    /// This is calulated the full Kalman Filter propagation of the track
    /// to the POCA and the projection of the propagated covarience matrix 
    /// along the line of the doca.
    inline float getDocaErr() const {
      return m_docaErr;
    }
    
    /// Return the POCA (in global coordinates)
    inline const Point& getPoca() const {
      return m_poca;
    }
    
    /// Return vector between the POCA and the closest edge of the detector element or GAP 
    inline const Vector& getPocaVector() const {
      return m_pocaVector;
    }
    
    
    /// set all the values
    inline void set(float arcLength, float doca, float docaErr, 
		    const Point& poca, const Vector& pocaVector);
    
    inline void set(const AcdPocaData& other);
    
    /// set individaul values
    inline void setArcLength(float val) { m_arcLength = val; }
    inline void setDoca(float val) { m_doca = val; }
    inline void setDocaErr(float val) { m_docaErr = val; }
    inline void setPoca(const Point& val) { m_poca = val; }
    inline void setPocaVector(const Vector& val) { m_pocaVector = val; }
    
    /// direct acces to members (to set them quickly)
    inline Point& accessPoca() { return m_poca; }
    inline Vector& accessPocaVector() { return m_pocaVector; }

    virtual void writeOut(MsgStream& stream ) const;
    
  protected:

    /// reset all the values to their default (lowercase to avoid conflict with TObject::Clear())
    void ini();
    
  private:
    
    /// The arclength at which the poca occurs
    float m_arcLength;
    
    /// The Distance of Closest Approach between track and tile or ribbon
    float m_doca;                 
    
    /// The Error on the DOCA, this is the track error ellipsoid project along the doca direction
    float m_docaErr;
    
    /// The Point of Closest Approach
    Point m_poca;
    
    /// The vector from the POCA to the Detector Element edge
    Vector m_pocaVector;
    
  };
    
};
#endif
