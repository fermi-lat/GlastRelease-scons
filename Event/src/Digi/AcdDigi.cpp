
//
// Inline code must be outside the class definition
//
#include "GlastEvent/Digi/AcdDigi.h"

using namespace idents;

/// Retrieve tile identifier
inline const ScintillatorId AcdDigi::ID() const
{
  return m_ID;
}

/// Retrieve pulse height
inline const int AcdDigi::PulseHeight() const
{
  return m_pulseHeight;
}

/// Retrieve veto
inline const bool AcdDigi::Veto() const
{
  return m_veto;
}

/// Retrieve low threshold
inline const bool AcdDigi::LowThreshold() const
{
  return m_lowThreshold;
}

/// Retrieve high threshold
inline const bool AcdDigi::HighThreshold() const
{
  return m_highThreshold;
}


