// $Header$

//
// Inline code must be outside the class definition
//
#include "math.h"
#include "GlastEvent/MonteCarlo/MCParticle.h"


/// Retrieve cell identifier
inline const VolumeID McPositionHit::volumeID() const
{
  return m_volumeID;
}
/// Update cell identifier
inline void McPositionHit::setVolumeID( VolumeID value )
{
  m_volumeID = value;
}


/// Retrieve entry member
inline const McPositionHit::HepPoint3D& entryPoint() const
{
  return m_entry;
}
inline McPositionHit::HepPoint3D& entryPoint()
{
  return m_entry;
}


/// Update Entry member
inline void McPositionHit::setEntryPoint( const HepPoint3D& value )
{
  m_entry = value;
}


/// Retrieve exit point
inline const HepPoint3D& McPositionHit::exitPoint() const
{
  return m_exit;
}
inline HepPoint3D& McPositionHit::exitPoint()
{
  return m_exit;
}


/// Update exit point
inline void McPositionHit::setExitPoint( const HepPoint3D& value )
{
  m_exit = value;
}


/// Retrieve deposited energy
inline double McPositionHit::depositedEnergy() const
{
  return m_depositedEnergy;
}


/// Update deposited energy
inline void McPositionHit::setDepositedEnergy( double value )
{
  m_depositedEnergy = value;
}

/// Retrieve depositing particle's energy
inline double McPositionHit::particleEnergy() const
{
  reuturn m_particleEnergy;
}

/// Set depositing particle's energy
inline void McPositionHit::setParticleEnergy(double value)
{
  m_particleEnergy = value;
}

/// Retrieve primary-origin flag
inline bool McPositionHit::primaryOrigin() const
{
  return m_packedFlags & ORIGIN_PRIMARY;
}
/// Update primary-origin flag
inline void McPositionHit::setPrimaryOrigin( bool value )
{
  if (value){
    m_packedFlags |= ORIGIN_PRIMARY;
  } else {
    m_packedFlags &= ~ORIGIN_PRIMARY;
  }
}
/// Retrieve calorimeter-shower-origin flag
inline bool McPositionHit::caloShowerOrigin() const
{
  return m_packedFlags & ORIGIN_CALOSHOWER;
}
/// Update calorimeter-shower-origin flag
inline void McPositionHit::setCaloShowerOrigin( bool value )
{
  if (value){
    m_packedFlags |= ORIGIN_CALOSHOWER;
  } else {
    m_packedFlags &= ~ORIGIN_CALOSHOWER;
  }
}

/// Retrieve whether this hit should be digitized
inline bool McPositionHit::needDigi() const
{
  return m_packedFlags & NEED_DIGI;
}
/// Update whether this hit should be digitized
inline void McPositionHit::setNeedDigi( bool value )
{
  if (value){
    m_packedFlags |= NEED_DIGI;
  } else {
    m_packedFlags &= ~NEED_DIGI;
  }
}

/// Retrieve hit's direction cosine
inline double McPositionHit::directionCosine() const
{
  double dx = m_exit.x()-m_entry.x();
  double dy = m_exit.y()-m_entry.y();
  double dz = m_exit.z()-m_entry.z();
  return dz / sqrt(dx * dx + dy * dy);
}

/// Retrieve member TOF
inline double McPositionHit::timeOfFlight() const
{
  return m_timeOfFlight;
}


/// Update TOF member
inline void McPositionHit::setTimeOfFlight( double value )
{
  m_timeOfFlight = value;
}


/// Retrieve pointer to MCParticle (const or non-const)
inline const MCParticle* McPositionHit::mcParticle() const
{
  return m_mcParticle; 
}
inline       MCParticle* McPositionHit::mcParticle()
{
  return m_mcParticle; 
}
/// Update pointer to MCParticle (by a C++ pointer or a smart reference)
inline void McPositionHit::setMCParticle( const MCParticle* value )
{
  m_mcParticle = value; 
  m_particleEnergy = value->fourMomentum().e();
}
inline void McPositionHit::setMCParticle( const SmartRef<MCParticle> value )
{ 
  m_mcParticle = value; 
  m_particleEnergy = value->fourMomentum().e();
}


/// Serialize the object for writing
inline StreamBuffer& McPositionHit::serialize( StreamBuffer& s ) const
{
  ContainedObject::serialize(s);
  return s
    << m_entry
    << m_depositedEnergy
    << m_timeOfFlight
    << m_mcParticle(this)
    << m_exit;
}


/// Serialize the object for reading
inline StreamBuffer& McPositionHit::serialize( StreamBuffer& s )
{
  ContainedObject::serialize(s);
  s >> m_entry
    >> m_depositedEnergy
    >> m_timeOfFlight
    >> m_mcParticle(this)
    >> m_exit;
  m_particleEnergy = m_mcParticle->fourMomentum().e();
  return s;
}


/// Fill the ASCII output stream
inline std::ostream& McPositionHit::fillStream( std::ostream& s ) const
{
  return s
    << "    base class McPositionHit :"
    << "\n        Entry point (x, y, z) = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_entry.x() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_entry.y() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_entry.z() << " )"
    << "\n        Deposited Energy      = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_depositedEnergy
    << "\n        Time of flight        = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_timeOfFlight
    << "\n        Exit point (x, y, z)  = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_exit.x() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_exit.y() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_exit.z() << " )";
}
