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
StreamBuffer& McParticle::serialize( StreamBuffer& s ) const
{
  ContainedObject::serialize(s);
  return s
    << m_particleID
    << m_particleProperty
    << m_subEvtID
    << m_statusFlags
    << m_originMcVertex(this)
    << m_endMcVertex(this);
}


/// Serialize the object for reading
StreamBuffer& McParticle::serialize( StreamBuffer& s )
{
  ContainedObject::serialize(s);
  s
    >> m_particleID
    >> m_particleProperty
    >> m_subEvtID
    >> m_statusFlags
    >> m_originMcVertex(this)
    >> m_endMcVertex(this);
  return s;
}


/// Fill the ASCII output stream
std::ostream& McParticle::fillStream( std::ostream& s ) const
{
  s << "class McParticle"
    << " (SubEvent:" << m_subEvtID << ")"
    << " :"
    << "\n    Particle ID                = " << m_particleID
    << "\n    Particle Property          = " << m_particleProperty
    << "\n    Sub Event ID               = " << m_subEvtID
    << "\n    Origin McVertex            = " << m_originMcVertex(this)
    << "\n    End McVertex               = " << m_endMcVertex(this);
  return s;
}
