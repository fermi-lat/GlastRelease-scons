// $Header$

#include <iostream>
#include <math.h>
#include "GlastEvent/MonteCarlo/McPositionHit.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"


// FIXME!!:
// The next `using' directive is an ad-hoc declaration for the transition
// phase to the namespace `GlastEvent'.
// After the transition is completed, it should be removed and each
// function definition should have the namespace identifier.
using namespace GlastEvent;

/// Retrieve hit's direction cosine
double McPositionHit::directionCosine() const
{
  double dx = m_exit.x()-m_entry.x();
  double dy = m_exit.y()-m_entry.y();
  double dz = m_exit.z()-m_entry.z();
  return dz / sqrt(dx * dx + dy * dy);
}


/// Serialize the object for writing
StreamBuffer& McPositionHit::serialize( StreamBuffer& s ) const
{
  ContainedObject::serialize(s);
  return s
    << m_volumeID
    << m_entry
    << m_exit
    << m_depositedEnergy
    << m_particleEnergy
    << m_timeOfFlight
    << m_mcParticle(this)
    << m_originMcParticle(this)
    << m_packedFlags;
}


/// Serialize the object for reading
StreamBuffer& McPositionHit::serialize( StreamBuffer& s )
{
  ContainedObject::serialize(s);
  return s
    >> m_volumeID
    >> m_entry
    >> m_exit
    >> m_depositedEnergy
    >> m_particleEnergy
    >> m_timeOfFlight
    >> m_mcParticle(this)
    >> m_originMcParticle(this)
    >> m_packedFlags;
}


/// Fill the ASCII output stream
std::ostream& McPositionHit::fillStream( std::ostream& s ) const
{
  return s
    << "    base class McPositionHit :"
    << "\n        Volume ID             = " << m_volumeID
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
    << "\n        Particle Energy       = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_particleEnergy
    << "\n        Time of flight        = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_timeOfFlight
    << "\n        Exit point (x, y, z)  = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_exit.x() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_exit.y() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_exit.z() << " )"
    << "\n        McParticle            = " << m_mcParticle(this)
    << "\n        ancestor McParticle   = " << m_originMcParticle(this);
}
