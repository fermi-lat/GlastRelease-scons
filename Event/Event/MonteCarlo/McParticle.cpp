// $Header$

//
// Inline code must be outside the class definition
//
#include "GlastEvent/MonteCarlo/MCVertex.h"



/// Retrieve particle identification
inline ParticleID McParticle::particleID() const
{
  return m_particleID;
}


/// Update particle identification
inline void McParticle::setParticleID( ParticleID value )
{
  m_particleID = value;
}


/// Retrieve particle property
inline ParticleProperty McParticle::particleProperty() const
{
  return m_particleProperty;
}


/// Update particle property
inline void McParticle::setParticleProperty( ParticleProperty value )
{
  m_particleProperty = value;
}


/// Retrieve whether this is a primary particle
inline bool McParticle::primaryParticle() const
{
  return m_statusFlags & STATUS_PRIMARY;
}


/// Set whether this is a primary particle
inline void McParticle::setPrimaryParticleFlag( bool value )
{
  if (value){
    m_statusFlags |= STATUS_PRIMARY;
  } else {
    m_statusFlags &= ~STATUS_PRIMARY;
  }
}


/// Retrieve pointer to origin vertex (const or non-const)
inline const MCVertex* McParticle::originMCVertex() const
{ 
  return m_originMCVertex; 
}
inline       MCVertex* McParticle::originMCVertex()
{ 
  return m_originMCVertex;
}


/// Update pointer to origin vertex (by a C++ pointer or a smart reference)
inline void McParticle::setOriginMCVertex( const MCVertex* value )
{ 
  m_originMCVertex = value; 
}
inline void McParticle::setOriginMCVertex( const SmartRef<MCVertex> value )
{ 
  m_originMCVertex = value; 
}


/// Retrieve pointer to vector of end vertex (const or non-const)
inline const MCVertex* McParticle::endMCVertex() const
{ 
  return m_endMCVertex; 
}
inline       MCVertex* McParticle::endMCVertex()
{ 
  return m_endMCVertex; 
}


/// Update pointer to end vertex (by a C++ pointer or a smart reference)
inline void McParticle::setEndMCVertex( const MCVertex* value )
{ 
  m_endMCVertex = value; 
}
inline void McParticle::setEndMCVertex( const SmartRef<MCVertex> value )
{ 
  m_endMCVertex = value; 
}


#ifdef NeedSubEvtID
/// Retrieve sub event ID 
inline short McParticle::subEvtID() const
{
  return m_subEvtID;
}


/// Set sub event ID 
inline void McParticle::setSubEvtID( short value )
{
  m_subEvtID = value;
}
#endif // NeedSubEvtID


/// Serialize the object for writing
inline StreamBuffer& McParticle::serialize( StreamBuffer& s ) const
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
inline StreamBuffer& McParticle::serialize( StreamBuffer& s )
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
inline std::ostream& McParticle::fillStream( std::ostream& s ) const
{
  s << "class McParticle"
#ifdef NeedSubEvtID
    << " (SubEvent:" << m_subEvtID << ")"
#endif // NeedSubEvtID
    << " :"
    << "\n    Particle ID                = " << m_particleID
    << "\n    Status Flags               = " << m_statusFlags;
  return s;
}
