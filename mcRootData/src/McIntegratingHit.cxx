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
    m_energyItem.clear();
    m_mcPartArr.Clear();
    m_energies.clear();
    m_packedFlags = 0;
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
    cout << "Size of EnergyMap: " << m_energyItem.size() << endl;
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


const McIntegratingHit::energyDepositMap& McIntegratingHit::getItemizedEnergy()
{
    // Purpose and Method:  Return a vector of pairs of McParticle* and energies.  
    //    We avoid storing a vector of pairs - due to the problems with storing 
    //    STL containers which contain a STL container.  Instead we store the map 
    //    as 2 separate arrays:  one contains the McParticles the other the 
    //    energies.  We construct the vector of pairs from these 2 arrays.  Note 
    //    that this vector of pairs is transient - it is never stored to a ROOT 
    //    file - it is constructed when needed.
    // Output:  An energyDepositMap - in this case a vector of pairs - where the
    //    pair consists of an McParticle* and a Double_t for the energy.
    if (m_energyItem.size() <= 0) {
        m_energyItem.clear();
        // TRefArray::Next does not seem to be implemented in 3.02 or 3.03.04
        // Just iterating over the entries in the arrays for now, rather than
        // using an iterator.
        //TRefArrayIter mcPartIter(&m_mcPartArr);
        McParticle *mcPart;
        UInt_t iEnergy = 0;
        //while( (mcPart = (McParticle*)mcPartIter.Next()) ) {
        for (iEnergy = 0; iEnergy < m_mcPartArr.GetEntries(); iEnergy++) {
            mcPart = (McParticle*) m_mcPartArr.At(iEnergy);
            m_energyItem.push_back( 
                std::pair<McParticle*, Double_t>(mcPart, m_energies[iEnergy]));
            //++iEnergy;
        }
    }

    return m_energyItem;
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

