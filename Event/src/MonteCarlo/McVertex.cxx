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
    << m_mcParticle(this)
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
    >> m_mcParticle(this)
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
    << "\n    Pair McParticle       = " << m_mcParticle(this)
    << "\n    Mother McParticle     = " << m_motherMcParticle(this)
    << "\n    Daughter McParticles  = ";
    SmartRefVector<McParticle>::const_iterator it;
    for (it = m_daughterMcParticles.begin(); it != m_daughterMcParticles.end(); it++){
      s << *it << "\n                            ";
    }
  return s;
}


/// Retrieve initial position
const HepPoint3D& McVertex::initialPosition () const
{
    return m_initialPosition;
}
/// Retrieve initial position
HepPoint3D& McVertex::initialPosition ()
{
    return m_initialPosition;
}
/// Update initial position
void McVertex::setInitialPosition (const HepPoint3D& value)
{
    m_initialPosition = value;
}

/// Retrieve final position
const HepPoint3D& McVertex::finalPosition () const
{
    return m_finalPosition;
}
/// Retrieve final position
HepPoint3D& McVertex::finalPosition ()
{
    return m_finalPosition;
}
/// Update final position
void McVertex::setFinalPosition (const HepPoint3D& value)
{
    m_finalPosition = value;
}

/// retrieve time of flight
double McVertex::timeOfFlight () const
{
    return m_timeOfFlight;
}
/// update time of flight
void McVertex::setTimeOfFlight (double value)
{
    m_timeOfFlight = value;
}
/// retrieve vertex type
McVertex::originType McVertex::vertexType () const
{
    return m_vertexType;
}
/// update vertex type
void McVertex::setVertexType (McVertex::originType value)
{
    m_vertexType = value;
}


/// Retrieve initial 4-momentum
const HepLorentzVector& McVertex::initialFourMomentum() const
{
    return m_initialFourMomentum;
}
/// Retrieve initial 4-momentum
HepLorentzVector& McVertex::initialFourMomentum()
{
    return m_initialFourMomentum;
}
/// Update initial 4-momentum
void McVertex::setInitialFourMomentum( const HepLorentzVector& value )
{
    m_initialFourMomentum = value;
}


/// Retrieve final 4-momentum
const HepLorentzVector& McVertex::finalFourMomentum() const
{
    return m_finalFourMomentum;
}
/// Retrieve final 4-momentum
HepLorentzVector& McVertex::finalFourMomentum()
{
    return m_finalFourMomentum;
}
/// Update final 4-momentum
void McVertex::setFinalFourMomentum( const HepLorentzVector& value )
{
    m_finalFourMomentum = value;
}


/// Retrieve pointer to the pair particle (const or non-const)
const McParticle* McVertex::mcParticle() const
{
    return m_mcParticle;
}
McParticle* McVertex::mcParticle()
{
    return m_mcParticle;
}
/// Update pointer to the pair particle (by a C++ pointer or a smart reference)
void McVertex::setMcParticle( McParticle* value )
{
    m_mcParticle = value;
}
void McVertex::setMcParticle( SmartRef<McParticle> value )
{
    m_mcParticle = value;
}


/// Retrieve pointer to mother particle (const or non-const)
const McParticle* McVertex::motherMcParticle() const
{
    return m_motherMcParticle;
}
McParticle* McVertex::motherMcParticle()
{
    return m_motherMcParticle;
}

/// Update pointer to mother particle (by a C++ pointer or a smart reference)
void McVertex::setMotherMcParticle( McParticle* value )
{
    m_motherMcParticle = value;
}
void McVertex::setMotherMcParticle( SmartRef<McParticle> value )
{
    m_motherMcParticle = value;
}

/// Retrieve pointer to vector of daughter particles (const or non-const)
const SmartRefVector<McParticle>& McVertex::daughterMcParticles() const
{
    return m_daughterMcParticles;
}
SmartRefVector<McParticle>& McVertex::daughterMcParticles()
{
    return m_daughterMcParticles;
}

/// Update all daughter particles
void McVertex::setDaughterMcParticles( const SmartRefVector<McParticle>& value )
{
    m_daughterMcParticles = value;
}

/// Remove all daughter particles
void McVertex::removeDaughterMcParticles()
{
    m_daughterMcParticles.clear();
}

/// Add single daughter particle to vector of daughter particles
///   (by a C++ pointer or a smart reference)
void McVertex::addDaughterMcParticle( McParticle* value )
{
    m_daughterMcParticles.push_back(value);
}
void McVertex::addDaughterMcParticle( SmartRef<McParticle> value )
{
    m_daughterMcParticles.push_back(value);
}


/// Retrieve sub event ID
short McVertex::subEvtID() const
{
    return m_subEvtID;
}
/// Set sub event ID
void McVertex::setSubEvtID( short value )
{
    m_subEvtID = value;
}
