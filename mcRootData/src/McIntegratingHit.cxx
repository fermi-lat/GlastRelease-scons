#include "mcRootData/McIntegratingHit.h"
#include <iostream>

ClassImp(McIntegratingHit)


McIntegratingHit::McIntegratingHit() :
m_packedFlags(0), m_totalEnergy(0.0)
{
    Clear();
}

McIntegratingHit::~McIntegratingHit() {
    
    Clear();    
}


void McIntegratingHit::initialize(const VolumeIdentifier& id) {
    m_volumeId = id;
}


void McIntegratingHit::Clear(Option_t *option)
{
    m_energyItem.clear();
    m_totalEnergy = 0.0;
    m_moment1Seed = TVector3(0., 0., 0.);
    m_moment2Seed = TVector3(0., 0., 0.);
    m_volumeId.Clear();
}

void McIntegratingHit::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    UInt_t p = 2;
    cout.precision(p);
    m_volumeId.Print(option);
    cout << "Flags: " << m_packedFlags 
        << " Energy: " << m_totalEnergy << endl;
    cout << "Mom1: (" << m_moment1Seed.X() << "," << m_moment1Seed.Y() << ","
        << m_moment1Seed.Z() << ")   ";
    cout << "Mom2: (" << m_moment2Seed.X() << "," << m_moment2Seed.Y() << ","
        << m_moment2Seed.Z() << ")" << endl;
}


const McIntegratingHit::energyDepositMap& McIntegratingHit::getItemizedEnergy() const
{
    return m_energyItem;
}


void McIntegratingHit::addEnergyItem(const Double_t& energy, McParticle* t, const TVector3& pos)
{
    m_energyItem[t] += energy; 
    TVector3 pos2 = TVector3(pos.X()*pos.X(), pos.Y()*pos.Y(), pos.Z()*pos.Z());
    m_totalEnergy += energy;
    m_moment1Seed = energy * pos;
    m_moment2Seed = energy * pos2;
}



const TVector3& McIntegratingHit::getMoment1 () const
{
    // Purpose and Method:  Retrieve the energy-weighted first moments of the
    //    position.
    return m_moment1Seed * (1./m_totalEnergy);
}


const TVector3& McIntegratingHit::getMoment2 () const
{
    // Purpose and Method:  Retrieve the energy-weighted second moments of the
    //    position
    return m_moment2Seed * (1./m_totalEnergy);
}

