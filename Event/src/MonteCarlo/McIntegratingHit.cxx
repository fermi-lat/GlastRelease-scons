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

namespace mc{
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
        const HepPoint3D& position  = it->first->finalPosition();
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


/// Retrieve volume identifier
const idents::VolumeIdentifier McIntegratingHit::volumeID() const
{
  return m_volumeID;
}


/// Update volume identifier
void McIntegratingHit::setVolumeID( idents::VolumeIdentifier value )
{
  m_volumeID = value;
}


/// Retrieve energy
double McIntegratingHit::totalEnergy() const
{
  return m_totalEnergy;
}


/// Retrieve the energy-weighted first moments of the position
const HepPoint3D McIntegratingHit::moment1 () const
{
    return m_moment1seed * (1./m_totalEnergy);
}
HepPoint3D McIntegratingHit::moment1 ()
{
    return m_moment1seed * (1./m_totalEnergy);
}


/// Retrieve the energy-weighted second moments of the position
const HepPoint3D McIntegratingHit::moment2 () const
{
    return m_moment2seed * (1./m_totalEnergy);
}
HepPoint3D McIntegratingHit::moment2 ()
{
    return m_moment2seed * (1./m_totalEnergy);
}


/// Retrieve itemized energy
const McIntegratingHit::energyDepositMap& McIntegratingHit::itemizedEnergy() const
{
  return m_energyItem;
}

McIntegratingHit::energyDepositMap& McIntegratingHit::itemizedEnergy()
{
  return m_energyItem;
}


/// Add an energyItem
void McIntegratingHit::addEnergyItem(const double& energy, mc::McParticle* t, const HepPoint3D& position)
{
    m_energyItem[t] += energy;

    HepPoint3D        position2 = HepPoint3D(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
    m_totalEnergy      += energy;
    m_moment1seed += energy * position;
    m_moment2seed += energy * position2;
}


/// Add an energyItem
void McIntegratingHit::addEnergyItem(const double& energy, SmartRef<mc::McParticle> t, const HepPoint3D& position)
{
    m_energyItem[t] += energy;

    HepPoint3D        position2 = HepPoint3D(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
    m_totalEnergy      += energy;
    m_moment1seed += energy * position;
    m_moment2seed += energy * position2;
}

}
