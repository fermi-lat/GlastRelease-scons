// $Header$

//
// Inline code must be outside the class definition
//
#include "GlastEvent/MonteCarlo/MCParticle.h"  // aka MCParticleKinematics


/// Retrieve position
inline const HepPoint3D& MCVertex::position () const
{
  return m_position;
}
inline HepPoint3D& MCVertex::position ()
{
  return m_position;
}
/// Update vertex position
inline void MCVertex::setPosition (const HepPoint3D& value)
{
  m_position = value;
}
/// retrieve time of flight
inline double MCVertex::timeOfFlight () const
{
  return m_timeOfFlight;
}
/// update time of flight
inline void MCVertex::setTimeOfFlight (double value)
{
  m_timeOfFlight = value;
}
/// retrieve vertex type
inline originType MCVertex::vertexType () const
{
  return t_vertexType;
}
/// update vertex type
inline void MCVertex::setVertexType (originType value)
{
  m_vertexType = value;
}


/// Retrieve pointer to mother particle (const or non-const)
inline const MCParticle* MCVertex::motherMCParticle() const                             {
  return m_motherMCParticle;
}
inline       MCParticle* MCVertex::motherMCParticle()
{
  return m_motherMCParticle;
}

/// Update pointer to mother particle (by a C++ pointer or a smart reference)
inline void MCVertex::setMotherMCParticle( const MCParticle* value )
{
  m_motherMCParticle = value;
}
inline void MCVertex::setMotherMCParticle( const SmartRef<MCParticle> value )
{
  m_motherMCParticle = value;
}

/// Retrieve pointer to vector of daughter particles (const or non-const)
inline const SmartRefVector<MCParticle>& MCVertex::daughterMCParticles() const
{
  return m_daughterMCParticles;
}
inline       SmartRefVector<MCParticle>& MCVertex::daughterMCParticles()
{
  return m_daughterMCParticles;
}

/// Update all daughter particles
inline void MCVertex::setDaughterMCParticles( const SmartRefVector<MCParticle>& value )
{
  m_daughterMCParticles = value;
}

/// Remove all daughter particles
inline void MCVertex::removeDaughterMCParticles()
{
  m_daughterMCParticles.clear();
}

/// Add single daughter particle to vector of daughter particles
///   (by a C++ pointer or a smart reference)
inline void MCVertex::addDaughterMCParticle( const MCParticle* value )
{
  m_daughterMCParticles.push_back(value);
}
inline void MCVertex::addDaughterMCParticle( const SmartRef<MCParticle> value )
{
  m_daughterMCParticles.push_back(value);
}


#ifdef NeedSubEvtID
/// Retrieve sub event ID
inline short MCVertex::subEvtID() const
{
  return m_subEvtID;
}
/// Set sub event ID
inline void MCVertex::setSubEvtID( short value )
{
  m_subEvtID = value;
}
#endif // NeedSubEvtID


/// Serialize the object for writing
inline StreamBuffer& MCVertex::serialize( StreamBuffer& s ) const
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
inline StreamBuffer& MCVertex::serialize( StreamBuffer& s )
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
inline std::ostream& MCVertex::fillStream( std::ostream& s ) const
{
  return s << "class MCVertex"
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
