#include "mcRootData/McIntegratingHit.h"

ClassImp(McIntegratingHit)


McIntegratingHit::McIntegratingHit() :
m_packedFlags(0), m_totalEnergy(0.0)
{
    Clear();
}

McIntegratingHit::~McIntegratingHit() {
    
    Clear();    
}


void McIntegratingHit::initialize(const VolumeIdentifier& id, Double_t e,
                            const TVector3 &moment1, const TVector3 &moment2) {
    
    //m_volumeId.initialize(id.getBits0to31(), id.getBits32to63(), id.size());
    m_volumeId = id;
    m_totalEnergy = e;
    m_moment1Seed.SetXYZ(moment1.X(), moment1.Y(), moment1.Z());
    m_moment2Seed.SetXYZ(moment2.X(), moment2.Y(), moment2.Z());
}


void McIntegratingHit::Clear(Option_t *option)
{
    m_energyItem.clear();
    m_totalEnergy = 0.0;
    m_moment1Seed = TVector3(0., 0., 0.);
    m_moment2Seed = TVector3(0., 0., 0.);
    m_volumeId.Clear();
}

const McIntegratingHit::energyDepositMap& McIntegratingHit::getItemizedEnergy() const
{
    return m_energyItem;
}


void McIntegratingHit::addEnergyItem(const Double_t& energy, McParticle* t)
{
    m_energyItem[t] += energy;    
}



const TVector3 McIntegratingHit::getMoment1 () const
{
    // Purpose and Method:  Retrieve the energy-weighted first moments of the
    //    position.
    return m_moment1Seed * (1./m_totalEnergy);
}


const TVector3 McIntegratingHit::getMoment2 () const
{
    // Purpose and Method:  Retrieve the energy-weighted second moments of the
    //    position
    return m_moment2Seed * (1./m_totalEnergy);
}

