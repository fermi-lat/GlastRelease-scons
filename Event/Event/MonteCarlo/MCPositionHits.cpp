// $Header$

//
// Inline code must be outside the class definition
//
#include "math.h"
#include "GlastEvent/MonteCarlo/MCParticle.h"


/// Retrieve cell identifier
inline const VolumeID MCPositionHits::volumeID() const
{
  return m_volumeID;
}
/// Update cell identifier
inline void MCPositionHits::setVolumeID( VolumeID value )
{
  m_volumeID = value;
}


/// Retrieve entry member
inline const MCPositionHits::HepPoint3D& entryPoint() const
{
  return m_entry;
}
inline MCPositionHits::HepPoint3D& entryPoint()
{
  return m_entry;
}


/// Update Entry member
inline void MCPositionHits::setEntryPoint( const HepPoint3D& value )
{
  m_entry = value;
}


/// Retrieve exit point
inline const HepPoint3D& MCPositionHits::exitPoint() const
{
  return m_exit;
}
inline HepPoint3D& MCPositionHits::exitPoint()
{
  return m_exit;
}


/// Update exit point
inline void MCPositionHits::setExitPoint( const HepPoint3D& value )
{
  m_exit = value;
}


/// Retrieve deposited energy
inline double MCPositionHits::depositedEnergy() const
{
  return m_depositedEnergy;
}


/// Update deposited energy
inline void MCPositionHits::setDepositedEnergy( double value )
{
  m_depositedEnergy = value;
}

/// Retrieve depositing particle's energy
inline double MCPositionHits::particleEnergy() const
{
  reuturn m_particleEnergy;
}

/// Set depositing particle's energy
inline void MCPositionHits::setParticleEnergy(double value)
{
  m_particleEnergy = value;
}

/// Retrieve primary-origin flag
inline bool MCPositionHits::primaryOrigin() const
{
  return m_packedFlags & ORIGIN_PRIMARY;
}
/// Update primary-origin flag
inline void MCPositionHits::setPrimaryOrigin( bool value )
{
  if (value){
    m_packedFlags |= ORIGIN_PRIMARY;
  } else {
    m_packedFlags &= ~ORIGIN_PRIMARY;
  }
}
/// Retrieve calorimeter-shower-origin flag
inline bool MCPositionHits::caloShowerOrigin() const
{
  return m_packedFlags & ORIGIN_CALOSHOWER;
}
/// Update calorimeter-shower-origin flag
inline void MCPositionHits::setCaloShowerOrigin( bool value )
{
  if (value){
    m_packedFlags |= ORIGIN_CALOSHOWER;
  } else {
    m_packedFlags &= ~ORIGIN_CALOSHOWER;
  }
}

/// Retrieve whether this hit should be digitized
inline bool MCPositionHits::needDigi() const
{
  return m_packedFlags & NEED_DIGI;
}
/// Update whether this hit should be digitized
inline void MCPositionHits::setNeedDigi( bool value )
{
  if (value){
    m_packedFlags |= NEED_DIGI;
  } else {
    m_packedFlags &= ~NEED_DIGI;
  }
}

/// Retrieve hit's direction cosine
inline double MCPositionHits::directionCosine() const
{
  double dx = m_exit.x()-m_entry.x();
  double dy = m_exit.y()-m_entry.y();
  double dz = m_exit.z()-m_entry.z();
  return dz / sqrt(dx * dx + dy * dy);
}

/// Retrieve member TOF
inline double MCPositionHits::timeOfFlight() const
{
  return m_timeOfFlight;
}


/// Update TOF member
inline void MCPositionHits::setTimeOfFlight( double value )
{
  m_timeOfFlight = value;
}


/// Retrieve pointer to MCParticle (const or non-const)
inline const MCParticle* MCPositionHits::mcParticle() const
{
  return m_mcParticle; 
}
inline       MCParticle* MCPositionHits::mcParticle()
{
  return m_mcParticle; 
}
/// Update pointer to MCParticle (by a C++ pointer or a smart reference)
inline void MCPositionHits::setMCParticle( const MCParticle* value )
{
  m_mcParticle = value; 
  m_particleEnergy = value->fourMomentum().e();
}
inline void MCPositionHits::setMCParticle( const SmartRef<MCParticle> value )
{ 
  m_mcParticle = value; 
  m_particleEnergy = value->fourMomentum().e();
}


/// Serialize the object for writing
inline StreamBuffer& MCPositionHits::serialize( StreamBuffer& s ) const
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
inline StreamBuffer& MCPositionHits::serialize( StreamBuffer& s )
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
inline std::ostream& MCPositionHits::fillStream( std::ostream& s ) const
{
  return s
    << "    base class MCPositionHits :"
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
