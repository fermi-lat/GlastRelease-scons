// $Header$

#include <iostream>
#include "CLHEP/Geometry/Point3D.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/Utilities/CLHEPStreams.h"


// FIXME!!:
// The next `using' directive is an ad-hoc declaration for the transition
// phase to the namespace `Event'.
// After the transition is completed, it should be removed and each
// function definition should have the namespace identifier.
//using namespace Event;

namespace Event{

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

void McIntegratingHit::setEnergyItems( double totEnergy, const double *energyArr,
                                      const HepPoint3D &moment1, const HepPoint3D &moment2) {
    // Purpose and Method:  Provide a mechanism to set the total energy and moments when
    //  reading from a persistent form, where we may not have the McParticles and their
    //  associated positions.

    m_totalEnergy = totEnergy;
    m_moment1seed = moment1 * m_totalEnergy;
    m_moment2seed = moment2 * m_totalEnergy;

    m_energyArray[0] = energyArr[0];
    m_energyArray[1] = energyArr[1];
    m_energyArray[2] = energyArr[2];
}

void McIntegratingHit::clearEnergyItems()
{
    m_energyItem.clear();
    m_energyItemId.clear();
    m_totalEnergy = 0.;
    m_moment1seed = HepPoint3D(0., 0., 0.);
    m_moment2seed = HepPoint3D(0., 0., 0.);
}


const idents::VolumeIdentifier McIntegratingHit::volumeID() const
{
  return m_volumeID;
}


void McIntegratingHit::setVolumeID( idents::VolumeIdentifier value )
{
  m_volumeID = value;
}


double McIntegratingHit::totalEnergy() const
{
  return m_totalEnergy;
}

const HepPoint3D McIntegratingHit::moment1 () const
{
    // Purpose and Method:   Retrieve the energy-weighted first moments of the
    //    position.
    return m_moment1seed * (1./m_totalEnergy);
}
HepPoint3D McIntegratingHit::moment1 ()
{
    return m_moment1seed * (1./m_totalEnergy);
}


const HepPoint3D McIntegratingHit::moment2 () const
{
    return m_moment2seed * (1./m_totalEnergy);
}

HepPoint3D McIntegratingHit::moment2 ()
{
    return m_moment2seed * (1./m_totalEnergy);
}


const McIntegratingHit::energyDepositMap& McIntegratingHit::itemizedEnergy() const
{
  return m_energyItem;
}

McIntegratingHit::energyDepositMap& McIntegratingHit::itemizedEnergy()
{
  return m_energyItem;
}

const McIntegratingHit::energyDepositMapId& McIntegratingHit::itemizedEnergyId() const
{
  return m_energyItemId;
}

McIntegratingHit::energyDepositMapId& McIntegratingHit::itemizedEnergyId()
{
  return m_energyItemId;
}

double McIntegratingHit::energyArray(Particle p) const { 
    return m_energyArray[p];
}

void McIntegratingHit::addEnergyItem(double energy, Particle p, const HepPoint3D& position)
{
    // Purpose and Method:  Add a McParticleId, energy pair to the collection.
    //    Update the total energy and moments.

    m_energyArray[p] += energy;

    HepPoint3D        position2 = HepPoint3D(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
    m_totalEnergy      += energy;
    m_moment1seed += energy * position;
    m_moment2seed += energy * position2;
}


void McIntegratingHit::addEnergyItem(const double& energy, Event::McParticle* t, const HepPoint3D& position)
{
    // Purpose and Method:  Add a McParticle*, energy pair to the collection.
    //    Update the total energy and moments.

    m_energyItem.push_back( std::pair<Event::McParticle*, double>(t, energy));

    HepPoint3D        position2 = HepPoint3D(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
    m_totalEnergy      += energy;
    m_moment1seed += energy * position;
    m_moment2seed += energy * position2;
}


void McIntegratingHit::addEnergyItem(const double& energy, SmartRef<Event::McParticle> t, const HepPoint3D& position)
{
    // Purpose and Method:  Add a McParticle*, energy pair to the collection.
    //    Update the total energy and moments.

    m_energyItem.push_back( std::pair<Event::McParticle*, double>(t, energy));
    HepPoint3D        position2 = HepPoint3D(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
    m_totalEnergy      += energy;
    m_moment1seed += energy * position;
    m_moment2seed += energy * position2;
}

}
