// $Header$

//
// Inline code must be outside the class definition
//
#include "GlastEvent/MonteCarlo/MCVertex.h"



/// Retrieve particle identification
inline ParticleID MCParticle::particleID() const
{
  return m_particleID;
}


/// Update particle identification
inline void MCParticle::setParticleID( ParticleID value )
{
  m_particleID = value;
}


/// Retrieve particle property
inline ParticleProperty MCParticle::particleProperty() const
{
  return m_particleProperty;
}


/// Update particle property
inline void MCParticle::setParticleProperty( ParticleProperty value )
{
  m_particleProperty = value;
}


/// Retrieve whether this is a primary particle
inline bool MCParticle::primaryParticle() const
{
  return m_statusFlags & STATUS_PRIMARY;
}


/// Set whether this is a primary particle
inline void MCParticle::setPrimaryParticleFlag( bool value )
{
  if (value){
    m_statusFlags |= STATUS_PRIMARY;
  } else {
    m_statusFlags &= ~STATUS_PRIMARY;
  }
}


/// Retrieve pointer to origin vertex (const or non-const)
inline const MCVertex* MCParticle::originMCVertex() const
{ 
  return m_originMCVertex; 
}
inline       MCVertex* MCParticle::originMCVertex()
{ 
  return m_originMCVertex;
}


/// Update pointer to origin vertex (by a C++ pointer or a smart reference)
inline void MCParticle::setOriginMCVertex( const MCVertex* value )
{ 
  m_originMCVertex = value; 
}
inline void MCParticle::setOriginMCVertex( const SmartRef<MCVertex> value )
{ 
  m_originMCVertex = value; 
}


/// Retrieve pointer to vector of end vertex (const or non-const)
inline const MCVertex* MCParticle::endMCVertex() const
{ 
  return m_endMCVertex; 
}
inline       MCVertex* MCParticle::endMCVertex()
{ 
  return m_endMCVertex; 
}


/// Update pointer to end vertex (by a C++ pointer or a smart reference)
inline void MCParticle::setEndMCVertex( const MCVertex* value )
{ 
  m_endMCVertex = value; 
}
inline void MCParticle::setEndMCVertex( const SmartRef<MCVertex> value )
{ 
  m_endMCVertex = value; 
}


#ifdef NeedSubEvtID
/// Retrieve sub event ID 
inline short MCParticle::subEvtID() const
{
  return m_subEvtID;
}


/// Set sub event ID 
inline void MCParticle::setSubEvtID( short value )
{
  m_subEvtID = value;
}
#endif // NeedSubEvtID


/// Serialize the object for writing
inline StreamBuffer& MCParticle::serialize( StreamBuffer& s ) const
{
  ContainedObject::serialize(s);
  return s
#ifdef NeedSubEvtID
    << m_subEvtID
#endif // NeedSubEvtID
    << m_particleID
    << m_statusFlags
    << m_originMCVertex(this)
    << m_endMCVertices(this);
}


/// Serialize the object for reading
inline StreamBuffer& MCParticle::serialize( StreamBuffer& s )
{
  ContainedObject::serialize(s);
  unsigned char u;
  s
#ifdef NeedSubEvtID
    >> m_subEvtID
#endif // NeedSubEvtID
    >> m_particleID
    >> m_statusFlags
    >> m_originMCVertex(this)
    >> m_endMCVertices(this);
  return s;
}


/// Fill the ASCII output stream
inline std::ostream& MCParticle::fillStream( std::ostream& s ) const
{
  s << "class MCParticle"
#ifdef NeedSubEvtID
    << " (SubEvent:" << m_subEvtID << ")"
#endif // NeedSubEvtID
    << " :"
    << "\n    Particle ID                = " << m_particleID
    << "\n    Status Flags               = " << m_statusFlags;
  return s;
}
