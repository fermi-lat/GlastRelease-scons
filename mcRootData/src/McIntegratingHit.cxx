#include "mcRootData/McIntegratingHit.h"
#include "TRefArray.h"
#include <iostream>

ClassImp(McIntegratingHit)


McIntegratingHit::McIntegratingHit() {
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
    m_mcPartArr.Clear();
    m_energies.clear();
    m_packedFlags = 0;
    m_totalEnergy = 0.0;
    m_moment1Seed = TVector3(0., 0., 0.);
    m_moment2Seed = TVector3(0., 0., 0.);
    m_volumeId.Clear();
    m_mapEntry = 0;
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
    cout << "Energy Map:  Size of energy vector: " << m_energies.size();
    cout << "    Size of McParticle array: " << m_mcPartArr.GetEntries() 
        << endl;
    TRefArrayIter mcPartIter(&m_mcPartArr);
    McParticle *mcPart = 0;
    UInt_t iEnergy = 0;
    cout << "Energy :" << endl;
    while( (mcPart = (McParticle*)mcPartIter.Next()) ) {
        cout << "( " << mcPart->getParticleId() << ", " << m_energies[iEnergy]
            << " )    " << endl;
        ++iEnergy;
    }

} 


const McParticle* McIntegratingHit::mapNext(Double_t &energy) {
    // Purpose and Method:  Access the next McParticle* and energy pair
    //   The transient data member, m_mapEntry, is used to access a specific
    //   pair.  If we have reached the end of the map, then a null McParticle*
    //   is returned.
    if ( m_mapEntry < m_mcPartArr.GetEntries() ) {
        energy = m_energies[m_mapEntry];
        return (McParticle*) m_mcPartArr.At(m_mapEntry++);
    }
    energy = 0;
    return 0;
}


void McIntegratingHit::addEnergyItem(const Double_t& energy, McParticle* t, 
                                     const TVector3& pos)
{
    // Purpose and Method:  Add a new McParticle, energy pair for a certain 
    //    position.  Stores the McParticle* and energy in the apporpriate 
    //    vectors.  Updates total energy and the first and second moments.
    m_mcPartArr.Add(t);
    m_energies.push_back(energy);
    TVector3 pos2 = TVector3(pos.X()*pos.X(), pos.Y()*pos.Y(), pos.Z()*pos.Z());
    m_totalEnergy += energy;
    m_moment1Seed = energy * pos;
    m_moment2Seed = energy * pos2;
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

