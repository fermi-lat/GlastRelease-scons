
//
// Inline code must be outside the class definition
//
#include "GlastEvent/Digi/AcdDigi.h"

using namespace idents;

/// Retrieve tile identifier
inline const ScintillatorId AcdDigi::getID() const
{
  return m_ID;
}

/// Retrieve pulse height
inline const int AcdDigi::getPulseHeight() const
{
  return m_pulseHeight;
}

/// Retrieve veto
inline const bool AcdDigi::getVeto() const
{
  return m_veto;
}

/// Retrieve low threshold
inline const bool AcdDigi::getLowThreshold() const
{
  return m_lowThreshold;
}

/// Retrieve high threshold
inline const bool AcdDigi::getHighThreshold() const
{
  return m_highThreshold;
}


/// Serialize the object for writing
inline StreamBuffer& AcdDigi::serialize( StreamBuffer& s ) const
{
  ContainedObject::serialize(s);
  return s
    << m_ID
    << m_packedPHA;
}


/// Serialize the object for reading
inline StreamBuffer& AcdDigi::serialize( StreamBuffer& s )
{
  ContainedObject::serialize(s);
 // s >> m_ID  need to add for ScintillatorId
   s >> m_packedPHA;

  // unpack the bits from m_packedPHA
  m_veto =  (m_packedPHA >> ADC_V_VETO) & ADC_M_VETO;
  m_lowThreshold =  (m_packedPHA >> ADC_V_LOWTHRESH) & ADC_M_LOWTHRESH;
  m_highThreshold =  (m_packedPHA >> ADC_V_HIGHTHRESH) & ADC_M_HIGHTHRESH;

  return s;
}


/// Fill the ASCII output stream
inline std::ostream& AcdDigi::fillStream( std::ostream& s ) const
{
  return s
    << "    base class AcdDigi :"
    << "\n        ID = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_ID << ", "
    << "\n        pulse height      = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_pulseHeight << ", "
    << "\n        veto              = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_veto   << " )"
    << "\n        Low Threshold     = "
    << m_lowThreshold << " )"
    << "\n        High Threshold     = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_highThreshold << " )"
    << "\n        packed PHA        = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_packedPHA << " )";
}

