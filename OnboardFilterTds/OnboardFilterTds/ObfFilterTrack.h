/** @file ObfFilterTrack.h
* @author Tracy Usher
*
* $Header$

*/
#ifndef ObfFilterTrack_H
#define ObfFilterTrack_H 1

#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/MsgStream.h"

static const CLID& CLID_ObfFilterTrack = InterfaceID("ObfFilterTrack", 1, 0);

namespace OnboardFilterTds {

/** 
* @class ObfFilterTrack
*
* @brief This TDS object contains the "best" track information from 
*        the onboard filter tracking.
*
*/
class ObfFilterTrack : virtual public DataObject
{
public:

    // Standard constructor
    ObfFilterTrack() : m_nXhits(0), m_nYhits(0), m_xInt(0), m_yInt(0), m_z(0), m_slpXZ(0), m_slpYZ(0) {}

    // Constructor with known parameters
    ObfFilterTrack(int xHits, int yHits, float xInt, float yInt, float z, float slpXZ, float slpYZ) :
        m_nXhits(xHits), m_nYhits(yHits), m_xInt(xInt), m_yInt(yInt), m_z(z), m_slpXZ(slpXZ), m_slpYZ(slpYZ) {}

    // TDS DataObject mumbo jumbo
    virtual const CLID& clID() const   { return ObfFilterTrack::classID(); }
    static const CLID& classID()       { return CLID_ObfFilterTrack; }

    // Standard destructor
    virtual ~ObfFilterTrack() {}

    // A special "fill" method for use in Root IO
    void initialize(int xHits, int yHits, float xInt, float yInt, float z, float slpXZ, float slpYZ);

    // Get the track paramters 
    const int   get_nXhits() const {return m_nXhits;}
    const int   get_nYhits() const {return m_nYhits;}
    const float get_xInt()   const {return m_xInt;}
    const float get_yInt()   const {return m_yInt;}
    const float get_z()      const {return m_z;}
    const float get_slpXZ()  const {return m_slpXZ;}
    const float get_slpYZ()  const {return m_slpYZ;}
        
    /// writes out the information of the cluster if msglevel is set to debug
    std::ostream& fillStream( std::ostream& s ) const;

private:
    int   m_nXhits;  // Number of hits on the X projection
    int   m_nYhits;  // Number of hits on the Y projection
    float m_xInt;    // Intercept of the X projection
    float m_yInt;    // Intercept of the Y projection
    float m_z;       // Z coordinate of intercept
    float m_slpXZ;   // slope X-Z plane
    float m_slpYZ;   // slope Y-Z plane
};

/// Fill the ASCII output stream
inline std::ostream& ObfFilterTrack::fillStream( std::ostream& s ) const
{
  return s << "class ObfFilterTrack"
    << " :"
    << "\n    Number hits X projection  = " << m_nXhits
    << "\n    Number hits Y projection  = " << m_nYhits
    << "\n    Intercept of X projection = " 
      << EventFloatFormat (EventFormat::width, EventFormat::precision ) 
      << m_xInt
    << "\n    Intercept of Y projection = " 
      << EventFloatFormat (EventFormat::width, EventFormat::precision ) 
      << m_yInt
    << "\n    Z Coordinate at intercept = " 
      << EventFloatFormat (EventFormat::width, EventFormat::precision ) 
      << m_z
    << "\n    Slope of X projection     = " 
      << EventFloatFormat (EventFormat::width, EventFormat::precision ) 
      << m_slpXZ
    << "\n    Slope of Y projection     = " 
      << EventFloatFormat (EventFormat::width, EventFormat::precision ) 
      << m_slpYZ
    << std::endl;
}

inline void ObfFilterTrack::initialize(int xHits, int yHits, float xInt, float yInt, float z, float slpXZ, float slpYZ)
{
    m_nXhits = xHits;
    m_nYhits = yHits;
    m_xInt   = xInt;
    m_yInt   = yInt;
    m_z      = z;
    m_slpXZ  = slpXZ;
    m_slpYZ  = slpYZ;
}


} // end namespace

#endif // ObfFilterTrack_H
