// $Header$

//
// Inline code must be outside the class definition
//
#include "GlastEvent/MonteCarlo/MCParticle.h"


/// Retrieve volume identifier
inline const VolumeID McIntegratingHit::volumeID() const
{
  return m_volumeID;
}


/// Update volume identifier
inline void McIntegratingHit::setVolumeID( VolumeID value )
{
  m_volumeID = value;
}


/// Retrieve energy
inline double McIntegratingHit::totalEnergy() const
{
  return m_totalEnergy;
}


/// Retrieve the energy-weighted first moments of the position
inline const HepPoint3D& McIntegratingHit::moment1 () const
{
    return m_moment1seed / m_totalEnergy;
}
inline HepPoint3D& McIntegratingHit::moment1 ()
{
    return m_moment1seed / m_totalEnergy;
}


/// Retrieve the energy-weighted second moments of the position
inline const HepPoint3D& McIntegratingHit::moment2 () const
{
    return m_moment2seed / m_totalEnergy;
}
inline HepPoint3D& McIntegratingHit::moment2 ()
{
    return m_moment2seed / m_totalEnergy;
}


/// Retrieve itemized energy
inline const energyDepositMap& McIntegratingHit::itemizedEnergy() const
{
  return m_energyItem;
}

inline       energyDepositMap& McIntegratingHit::itemizedEnergy()
{
  return m_energyItem;
}


/// Update all energyInfos
inline void McIntegratingHit::setEnergyItems( const energyDepositMap& value )
{
    m_energyItem = value;
    m_totalEnergy = 0.;
    m_moment1seed = HepPoint3D(0., 0., 0.);
    m_moment2seed = HepPoint3D(0., 0., 0.);
    typedef energyDepositMap::const_iterator CI;
    for (CI it = m_energyItem.begin(); it != m_energyItem.end(); it++){
        double&      energy    = it->first;
        HepPoint3D&  position  = it->second->endMCVertex()->position();
        HepPoint3D   position2 = HepPoint3D(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
        m_totalEnergy += energy;
        m_moment1seed += energy * position;
        m_moment2seed += energy * position2;
    }
}


/// Remove all energyInfos
inline void McIntegratingHit::clearEnergyItems()
{
    m_energyItem.clear();
    m_totalEnergy = 0.;
    m_moment1seed = HepPoint3D(0., 0., 0.);
    m_moment2seed = HepPoint3D(0., 0., 0.);
}


/// Add an energyItem
inline void McIntegratingHit::addEnergyItem( const double E, const MCParticle* t)
{
    m_energyItem[t] = E;

    double&      energy    = E;
    HepPoint3D&  position  = t->endMCVertex()->position();
    HepPoint3D   position2 = HepPoint3D(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
    m_totalEnergy      += energy;
    m_moment1seed += energy * position;
    m_moment2seed += energy * position2;
}


/// Add an energyItem
inline void McIntegratingHit::addEnergyItem( const double E, const SmartRef<MCParticle>& t)
{
    m_energyItem[t] = E;

    double&      energy    = E;
    HepPoint3D&  position  = t->endMCVertex()->position();
    HepPoint3D   position2 = HepPoint3D(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
    m_totalEnergy      += energy;
    m_moment1seed += energy * position;
    m_moment2seed += energy * position2;
}


/// Retrieve primary-origin flag
inline bool McIntegratingHit::primaryOrigin() const
{
    return m_packedFlags & ORIGIN_PRIMARY;
}


/// Update primary-origin flag
inline void McIntegratingHit::setPrimaryOrigin( bool value )
{
    if (value){
        m_packedFlags |= ORIGIN_PRIMARY;
    } else {
        m_packedFlags &= ~ORIGIN_PRIMARY;
    }
}


/// Retrieve whether this hit should be digitized
inline bool McIntegratingHit::primaryOrigin() const
{
    return m_packedFlags & NEED_DIGI;
}


/// Update whether this hit should be digitized
inline void McIntegratingHit::setPrimaryOrigin( bool value )
{
    if (value){
        m_packedFlags |= NEED_DIGI;
    } else {
        m_packedFlags &= ~NEED_DIGI;
    }
}


/// Serialize the object for writing
inline StreamBuffer& McIntegratingHit::serialize( StreamBuffer& s ) const
{
    ContainedObject::serialize(s);
    return s
      << m_volumeID
      << m_energyItem(this)
      << m_packedFlags;
}


/// Serialize the object for reading
inline StreamBuffer& McIntegratingHit::serialize( StreamBuffer& s )
{
    ContainedObject::serialize(s);
    s >> m_volumeID
      >> m_energyItem(this)
      >> m_packedFlags;

    m_totalEnergy = 0.;
    m_moment1seed = HepPoint3D(0., 0., 0.);
    m_moment2seed = HepPoint3D(0., 0., 0.);
    typedef energyDepositMap::const_iterator CI;
    for (CI it = m_energyItem.begin(); it != m_energyItem.end(); it++){
        double&      energy    = it->first;
        HepPoint3D&  position  = it->second->endMCVertex()->position();
        HepPoint3D   position2 = HepPoint3D(position.x()*position.x(), position.
y()*position.y(), position.z()*position.z());
        m_totalEnergy += energy;
        m_moment1seed += energy * position;
        m_moment2seed += energy * position2;
    }

    return s;
}


/// Fill the ASCII output stream
inline std::ostream& McIntegratingHit::fillStream( std::ostream& s ) const
{
    s << "class MCCaloHitBase :"
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
