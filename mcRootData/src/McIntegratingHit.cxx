#include "mcRootData/McIntegratingHit.h"

ClassImp(McIntegratingHit)


McIntegratingHit::McIntegratingHit()
{
    
}

McIntegratingHit::~McIntegratingHit() {
    
    Clear();
    
}


void McIntegratingHit::initialize(VolumeIdentifier id, Double_t e,
                            const TVector3 &moment1, const TVector3 &moment2) {
    
    m_volumeId = id;
    m_totalEnergy = e;
    m_moment1Seed = moment1;
    m_moment2Seed = moment2;
}


void McIntegratingHit::Clear(Option_t *option)
{
    clearEnergyItems();
}

const McIntegratingHit::energyDepositMap& McIntegratingHit::getItemizedEnergy() const
{
    return m_energyItem;
}


void McIntegratingHit::addEnergyItem(const Double_t& energy, McParticle* t, const TVector3& position)
{
    m_energyItem[t] += energy;
    
    TVector3 position2 = TVector3(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
    m_totalEnergy += energy;
    m_moment1Seed += energy * position;
    m_moment2Seed += energy * position2;
}

void McIntegratingHit::setEnergyItems( const energyDepositMap& value )
{
    m_energyItem = value;
    m_totalEnergy = 0.;
    m_moment1Seed = TVector3(0., 0., 0.);
    m_moment2Seed = TVector3(0., 0., 0.);
    typedef energyDepositMap::const_iterator CI;
    for (CI it = m_energyItem.begin(); it != m_energyItem.end(); it++){
        const Double_t&     energy    = it->second;
        const TVector3& position  = it->first->getFinalPosition();
        TVector3 position2 = TVector3(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
        m_totalEnergy += energy;
        m_moment1Seed += energy * position;
        m_moment2Seed += energy * position2;
    }
}


/// Remove all energyInfos
void McIntegratingHit::clearEnergyItems()
{
    m_energyItem.clear();
    m_totalEnergy = 0.;
    m_moment1Seed = TVector3(0., 0., 0.);
    m_moment2Seed = TVector3(0., 0., 0.);
}


/// Retrieve the energy-weighted first moments of the position
const TVector3 McIntegratingHit::getMoment1 () const
{
    return m_moment1Seed * (1./m_totalEnergy);
}


/// Retrieve the energy-weighted second moments of the position
const TVector3 McIntegratingHit::getMoment2 () const
{
    return m_moment2Seed * (1./m_totalEnergy);
}


