// $Header$

#include <iostream>
#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"


// FIXME!!:
// The next `using' directive is an ad-hoc declaration for the transition
// phase to the namespace `GlastEvent'.
// After the transition is completed, it should be removed and each
// function definition should have the namespace identifier.
using namespace GlastEvent;

/// Serialize the object for writing

/// Retrieve particle identification
ParticleID McParticle::particleID() const
{
  return m_particleID;
}


/// Update particle identification
void McParticle::setParticleID( ParticleID value )
{
  m_particleID = value;
}


/// Retrieve particle property
McParticle::StdHepId McParticle::particleProperty() const
{
  return m_particleProperty;
}


/// Update particle property
void McParticle::setParticleProperty( McParticle::StdHepId value )
{
  m_particleProperty = value;
}


/// Retrieve whether this is a primary particle
bool McParticle::primaryParticle() const
{
  using GlastEvent::McConstants::PRIMARY;
  return m_statusFlags & PRIMARY;
}


/// Set whether this is a primary particle
void McParticle::setPrimaryParticleFlag( bool value )
{
  using GlastEvent::McConstants::PRIMARY;
  if (value){
    m_statusFlags |= PRIMARY;
  } else {
    m_statusFlags &= ~PRIMARY;
  }
}


/// Retrieve pointer to the vertex (const or non-const)
const McVertex* McParticle::mcVertex() const
{ 
  return m_mcVertex; 
}
      McVertex* McParticle::mcVertex()
{ 
  return m_mcVertex;
}


/// Update pointer to the vertex (by a C++ pointer or a smart reference)
void McParticle::setMcVertex( McVertex* value )
{ 
  m_mcVertex = value; 
}
void McParticle::setMcVertex( SmartRef<McVertex> value )
{ 
  m_mcVertex = value; 
}


/// Retrieve sub event ID 
short McParticle::subEvtID() const
{
  return m_subEvtID;
}


/// Set sub event ID 
void McParticle::setSubEvtID( short value )
{
  m_subEvtID = value;
}
