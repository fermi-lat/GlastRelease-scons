// $Header$

//
// Inline code must be outside the class definition
//
#include "GlastEvent/MonteCarlo/MCParticle.h"  // aka MCParticleKinematics


/// Retrieve initial position
inline const HepPoint3D& McVertex::initialPosition () const
{
  return m_initialPosition;
}
/// Retrieve initial position
inline HepPoint3D& McVertex::initialPosition ()
{
  return m_initialPosition;
}
/// Update initial position
inline void McVertex::setInitialPosition (const HepPoint3D& value)
{
  m_initialPosition = value;
}

/// Retrieve final position
inline const HepPoint3D& McVertex::finalPosition () const
{
  return m_finalPosition;
}
/// Retrieve final position
inline HepPoint3D& McVertex::finalPosition ()
{
  return m_finalPosition;
}
/// Update final position
inline void McVertex::setFinalPosition (const HepPoint3D& value)
{
  m_finalPosition = value;
}

/// retrieve time of flight
inline double McVertex::timeOfFlight () const
{
  return m_timeOfFlight;
}
/// update time of flight
inline void McVertex::setTimeOfFlight (double value)
{
  m_timeOfFlight = value;
}
/// retrieve vertex type
inline originType McVertex::vertexType () const
{
  return t_vertexType;
}
/// update vertex type
inline void McVertex::setVertexType (originType value)
{
  m_vertexType = value;
}


/// Retrieve pointer to mother particle (const or non-const)
inline const MCParticle* McVertex::motherMCParticle() const                             {
  return m_motherMCParticle;
}
inline       MCParticle* McVertex::motherMCParticle()
{
  return m_motherMCParticle;
}

/// Update pointer to mother particle (by a C++ pointer or a smart reference)
inline void McVertex::setMotherMCParticle( const MCParticle* value )
{
  m_motherMCParticle = value;
}
inline void McVertex::setMotherMCParticle( const SmartRef<MCParticle> value )
{
  m_motherMCParticle = value;
}

/// Retrieve pointer to vector of daughter particles (const or non-const)
inline const SmartRefVector<MCParticle>& McVertex::daughterMCParticles() const
{
  return m_daughterMCParticles;
}
inline       SmartRefVector<MCParticle>& McVertex::daughterMCParticles()
{
  return m_daughterMCParticles;
}

/// Update all daughter particles
inline void McVertex::setDaughterMCParticles( const SmartRefVector<MCParticle>& value )
{
  m_daughterMCParticles = value;
}

/// Remove all daughter particles
inline void McVertex::removeDaughterMCParticles()
{
  m_daughterMCParticles.clear();
}

/// Add single daughter particle to vector of daughter particles
///   (by a C++ pointer or a smart reference)
inline void McVertex::addDaughterMCParticle( const MCParticle* value )
{
  m_daughterMCParticles.push_back(value);
}
inline void McVertex::addDaughterMCParticle( const SmartRef<MCParticle> value )
{
  m_daughterMCParticles.push_back(value);
}


#ifdef NeedSubEvtID
/// Retrieve sub event ID
inline short McVertex::subEvtID() const
{
  return m_subEvtID;
}
/// Set sub event ID
inline void McVertex::setSubEvtID( short value )
{
  m_subEvtID = value;
}
#endif // NeedSubEvtID


/// Serialize the object for writing
inline StreamBuffer& McVertex::serialize( StreamBuffer& s ) const
{
  ContainedObject::serialize(s);
  return s
#ifdef NeedSubEvtID
    << m_subEvtID
#endif // NeedSubEvtID
    << m_position
    << m_timeOfFlight
    << m_vertexType
    << m_initialFourMomentum
    << m_finalFourMomentum
    << m_motherMCParticle(this)
    << m_daughterMCParticles(this);
}

/// Serialize the object for reading
inline StreamBuffer& McVertex::serialize( StreamBuffer& s )
{
  int tmp;
  ContainedObject::serialize(s);
  s
#ifdef NeedSubEvtID
    >> m_subEvtID
#endif // NeedSubEvtID
    >> m_position
    >> m_timeOfFlight
    >> m_vertexType
    >> m_initialFourMomentum
    >> m_finalFourMomentum
    >> m_motherMCParticle(this)
    >> m_daughterMCParticles(this);

  return s;
}

/// Fill the ASCII output stream
inline std::ostream& McVertex::fillStream( std::ostream& s ) const
{
  return s << "class McVertex"
#ifdef NeedSubEvtID
    << "  (SubEvent:" << m_subEvtID << ")"
#endif // NeedSubEvtID
    << " :"
    << "\n    Position (x, y, z) = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_position.x() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_position.y() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_position.z() << " )"
    << "\n    Time of flight     = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_timeOfFlight
    << "\n    Initial 4-momentum (px, py, pz, E) = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_initialFourMomentum.px() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_initialFourMomentum.py() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_initialFourMomentum.pz() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_initialFourMomentum.e()  << " )"
    << "\n    Final 4-momentum (px, py, pz, E) = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_finalFourMomentum.px() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_finalFourMomentum.py() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_finalFourMomentum.pz() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_finalFourMomentum.e()  << " )"
    << "\n    Vertex Type        = "
    << int(m_vertexType);
}
