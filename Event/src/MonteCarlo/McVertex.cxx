// $Header$

#include <iostream>
#include "GlastEvent/MonteCarlo/McVertex.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"


// FIXME!!:
// The next `using' directive is an ad-hoc declaration for the transition
// phase to the namespace `GlastEvent'.
// After the transition is completed, it should be removed and each
// function definition should have the namespace identifier.
using namespace GlastEvent;

/// Serialize the object for writing
StreamBuffer& McVertex::serialize( StreamBuffer& s ) const
{
  int tmp = m_vertexType;
  ContainedObject::serialize(s);
  return s
    << m_subEvtID
    << m_initialPosition
    << m_finalPosition
    << m_timeOfFlight
    << tmp // m_vertexType
    << m_initialFourMomentum
    << m_finalFourMomentum
    << m_motherMcParticle(this)
    << m_daughterMcParticles(this);
}

/// Serialize the object for reading
StreamBuffer& McVertex::serialize( StreamBuffer& s )
{
  int tmp;
  ContainedObject::serialize(s);
  s
    >> m_subEvtID
    >> m_initialPosition
    >> m_finalPosition
    >> m_timeOfFlight
    >> tmp // m_vertexType
    >> m_initialFourMomentum
    >> m_finalFourMomentum
    >> m_motherMcParticle(this)
    >> m_daughterMcParticles(this);

  m_vertexType = originType(tmp);
  return s;
}

/// Fill the ASCII output stream
std::ostream& McVertex::fillStream( std::ostream& s ) const
{
  s << "class McVertex"
    << "  (SubEvent:" << m_subEvtID << ")"
    << " :"
    << "\n    initialPosition (x, y, z) = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_initialPosition.x() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_initialPosition.y() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_initialPosition.z() << " )"
    << "\n    finalPosition (x, y, z) = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_finalPosition.x() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_finalPosition.y() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_finalPosition.z() << " )"
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
    << "\n    Vertex Type           = " << int(m_vertexType)
    << "\n    Mother McParticle     = " << m_motherMcParticle(this)
    << "\n    Daughter McParticles  = ";
    SmartRefVector<McParticle>::const_iterator it;
    for (it = m_daughterMcParticles.begin(); it != m_daughterMcParticles.end(); it++){
      s << *it << "\n                            ";
    }
  return s;
}
