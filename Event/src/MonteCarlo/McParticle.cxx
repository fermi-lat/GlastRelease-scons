// $Header$

#include <iostream>
#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"


namespace mc {

/// Retrieve particle property
McParticle::StdHepId McParticle::particleProperty() const
{
  return m_particleID;
}


/// Update particle property
void McParticle::setParticleProperty( McParticle::StdHepId value )
{
  m_particleID = value;
}


/// Retrieve whether this is a primary particle
bool McParticle::primaryParticle() const
{
  return m_statusFlags & PRIMARY;
}




} // namespace



