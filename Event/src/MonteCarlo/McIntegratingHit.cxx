// $Header$

#include <iostream>
#include "CLHEP/Geometry/Point3D.h"
#include "GlastEvent/MonteCarlo/McIntegratingHit.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"


// FIXME!!:
// The next `using' directive is an ad-hoc declaration for the transition
// phase to the namespace `GlastEvent'.
// After the transition is completed, it should be removed and each
// function definition should have the namespace identifier.
using namespace GlastEvent;

/// Update all energyInfos
void McIntegratingHit::setEnergyItems( const energyDepositMap& value )
{
    m_energyItem = value;
    m_totalEnergy = 0.;
    m_moment1seed = HepPoint3D(0., 0., 0.);
    m_moment2seed = HepPoint3D(0., 0., 0.);
    typedef energyDepositMap::const_iterator CI;
    for (CI it = m_energyItem.begin(); it != m_energyItem.end(); it++){
        const double&     energy    = it->second;
        const HepPoint3D& position  = it->first->mcVertex()->finalPosition();
        HepPoint3D        position2 = HepPoint3D(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
        m_totalEnergy += energy;
        m_moment1seed += energy * position;
        m_moment2seed += energy * position2;
    }
}


/// Remove all energyInfos
void McIntegratingHit::clearEnergyItems()
{
    m_energyItem.clear();
    m_totalEnergy = 0.;
    m_moment1seed = HepPoint3D(0., 0., 0.);
    m_moment2seed = HepPoint3D(0., 0., 0.);
}


/// Serialize the object for writing
StreamBuffer& McIntegratingHit::serialize( StreamBuffer& s ) const
{
    ContainedObject::serialize(s);
    s
      << m_volumeID
      << m_totalEnergy
      << m_moment1seed
      << m_moment2seed
//    << m_energyItem(this)	// The operator<< has not implemented yet. FIX ME!!
      << m_energyItem.size();
    energyDepositMap::const_iterator it;
    for (it = m_energyItem.begin(); it != m_energyItem.end(); it++){
        s << it->first(this)
          << it->second;
    }
    return s
      << m_packedFlags;
}


/// Serialize the object for reading
StreamBuffer& McIntegratingHit::serialize( StreamBuffer& s )
{
    energyDepositMap::size_type m_energyItem_size;
    ContainedObject::serialize(s);
    s
      >> m_volumeID
      >> m_totalEnergy
      >> m_moment1seed
      >> m_moment2seed
//    >> m_energyItem(this)	// The operator<< has not implemented yet. FIX ME!!
      >> m_energyItem_size;
    m_energyItem.clear();
    energyDepositMap::size_type i;
    for (i = 0; i < m_energyItem_size; i++){
        SmartRef<McParticle> first;
        double               second;
        s >> first(this)
          >> second;
        m_energyItem[first] = second;
    }
    return s
      >> m_packedFlags;
}


/// Fill the ASCII output stream
std::ostream& McIntegratingHit::fillStream( std::ostream& s ) const
{
    s << "class McCaloHitBase :"
      << "\n    Deposited Energy        = "
      << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
      << m_totalEnergy
      << "\n    First moment (x, y, z)  = "
      << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
      << m_moment1seed.x() / m_totalEnergy << ", "
      << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
      << m_moment1seed.y() / m_totalEnergy << ", "
      << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
      << m_moment1seed.z() / m_totalEnergy << " )"
      << "\n    Second moment (x, y, z) = "
      << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
      << m_moment2seed.x() / m_totalEnergy << ", "
      << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
      << m_moment2seed.y() / m_totalEnergy << ", "
      << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
      << m_moment2seed.z() / m_totalEnergy << " )"
      << "\n    Volume ID               = " << m_volumeID;
    return s;
}
