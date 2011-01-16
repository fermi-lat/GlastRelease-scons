#include <mcRootData/McIntegratingHit.h>
#include <commonRootData/RootDataUtil.h>
#include "TRefArray.h"
#include "TObjString.h"
#include "Riostream.h"

#include "McObjectManager.h"

ClassImp(McIntegratingHit)
ClassImp(McXtalEnergyDep)


McIntegratingHit::McIntegratingHit() {
    Clear();
    m_xtalEnergyMap.SetOwnerKeyValue();
    m_xtalEnergyMap.SetOwnerValue();
}

McIntegratingHit::~McIntegratingHit() {
    
    Clear();    
}


void McIntegratingHit::initialize(const VolumeIdentifier& id) {
    m_volumeId = id;
}


void* McIntegratingHit::operator new(size_t /*size*/)
{
    McIntegratingHit* temp = McObjectManager::getPointer()->getNewMcIntegratingHit();

    // Since we are re-using, make sure some of the data members which might allocate memory are cleaned up
    temp->m_mcPartArr.Clear();
    temp->m_energyPtrArr.clear();

    return temp;
}

void* McIntegratingHit::operator new(size_t /*size*/, void* vp)
{
    return vp;
}

McIntegratingHit& McIntegratingHit::operator =(const McIntegratingHit& rhs)
{
    m_totalEnergy    = rhs.m_totalEnergy;
    m_energyArray[0] = rhs.m_energyArray[0];
    m_energyArray[1] = rhs.m_energyArray[1];
    m_energyArray[2] = rhs.m_energyArray[2];
    m_packedFlags    = rhs.m_packedFlags;
    m_volumeId       = rhs.m_volumeId;
    m_moment1Seed    = rhs.m_moment1Seed;
    m_moment2Seed    = rhs.m_moment2Seed;
    m_mapPtr         = rhs.m_mapPtr; 

    for(int idx=0; idx < rhs.m_mcPartArr.GetEntries(); idx++)
        m_mcPartArr.Add(rhs.m_mcPartArr.At(idx));

    return *this;
}

void McIntegratingHit::Clear( Option_t * )
{
    m_mcPartArr.Clear();
    m_energyPtrArr.clear();
    m_packedFlags = 0;
    m_totalEnergy = 0.0;
    m_moment1Seed = TVector3(0., 0., 0.);
    m_moment2Seed = TVector3(0., 0., 0.);
    m_volumeId.Clear();
    m_mapPtr = 0;
    m_energyArray[0] = 0.; m_energyArray[1] = 0.; m_energyArray[2] = 0.;
    m_xtalEnergyMap.Clear();
}

// dummy data, just for tests
void McIntegratingHit::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

    Clear() ;
    
    VolumeIdentifier id ;
    id.Clear() ;
    id.append(0) ;
    initialize(id) ;
    
    //TVector3 pos = mcPart->getFinalPosition();
    //TVector3 pos2(1.3, 0.0, 12.0);
    //addEnergyItem(1.5, mcPart, pos);
    
    double totE = 5.5;
    double energyArr[3] = { 2.5, 3.0, 0.0 };
    TVector3 moment1(1.0, 2.0, 3.0);
    TVector3 moment2(2.0, 4.0, 6.0);
    setEnergyItems(totE,energyArr,moment1,moment2) ;

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att,ref.get ## att,#att)

Bool_t McIntegratingHit::CompareInRange( const McIntegratingHit & ref, const std::string & name ) const {

    Bool_t result = true ;
    
    result = COMPARE_IN_RANGE(VolumeId()) && result ;
    
    result = COMPARE_IN_RANGE(TotalEnergy()) && result ;
    result = COMPARE_IN_RANGE(Moment1()) && result ;
    result = COMPARE_IN_RANGE(Moment2()) && result ;

    result = COMPARE_IN_RANGE(McParticleEnergy(McIntegratingHit::PRIMARY)) && result ;
    result = COMPARE_IN_RANGE(McParticleEnergy(McIntegratingHit::ELECTRON)) && result ;
    result = COMPARE_IN_RANGE(McParticleEnergy(McIntegratingHit::POSITRON)) && result ;

    itemizedEnergyReset() ;
    ref.itemizedEnergyReset() ;
    result = rootdatautil::CompareInRange(itemizedEnergySize(),ref.itemizedEnergySize(),"itemizedEnergySize") && result ;  
    const McParticle * myPart, * refPart ;
    double myEnergy, refEnergy ;
    int count = 0 ;
    while ( (myPart = itemizedEnergyNext(myEnergy)) ) {
      ++count ;
      refPart = ref.itemizedEnergyNext(refEnergy) ;
      result = rootdatautil::CompareInRange(*myPart,*refPart,"itemizedEnergyNext Particle") && result ;
      result = rootdatautil::CompareInRange(myEnergy,refEnergy,"itemizedEnergyNext Energy") && result ;
    }
    if (count==0) {
      std::cout
        <<"No Error : cannot read map which stores TRefs using " 
        <<" compiled code and ROOT 3.02.07 "
        <<std::endl ;
    }
      
    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

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
    if (m_energyPtrArr.size() > 0) {
        cout << "Energy Map:  Size of energy vector: " << m_energyPtrArr.size();
        cout << "    Size of McParticle array: " << m_mcPartArr.GetEntries() 
        << endl;
    } else {
        cout << "energies stored: " << m_energyArray[0] << " , "
            << m_energyArray[1] << " , " << m_energyArray[2] << endl;
    }
    TRefArrayIter mcPartIter(&m_mcPartArr);
    McParticle *mcPart = 0;
    UInt_t iEnergy = 0;
    cout << "Energy :" << endl;
    while( (mcPart = (McParticle*)mcPartIter.Next()) ) {
        cout << "( " << mcPart->getParticleId() << ", " << m_energyPtrArr[iEnergy]
            << " )    " << endl;
        ++iEnergy;
    }

} 


const McParticle* McIntegratingHit::itemizedEnergyNext(Double_t &energy) const {
    // Purpose and Method:  Access the next McParticle* and energy pair
    //   The transient data member, m_mapEntry, is used to access a specific
    //   pair.  If we have reached the end of the map, then a null McParticle*
    //   is returned.
    if ( m_mapPtr < (UInt_t)m_mcPartArr.GetEntries() ) {
        energy = m_energyPtrArr[m_mapPtr];
        return (McParticle*) m_mcPartArr.At(m_mapPtr++);
    }
    energy = 0;
    return 0;
}

Double_t McIntegratingHit::getMcParticleEnergy(Particle p) const {
    return m_energyArray[p];
}

void McIntegratingHit::addEnergyItem(const Double_t& energy, McParticle* t, 
                                     const TVector3& pos)
{
    // Purpose and Method:  Add a new McParticle, energy pair for a certain 
    //    position.  Stores the McParticle* and energy in the apporpriate 
    //    vectors.  Updates total energy and the first and second moments.
    m_mcPartArr.Add(t);
    m_energyPtrArr.push_back(energy);
    TVector3 pos2 = TVector3(pos.X()*pos.X(), pos.Y()*pos.Y(), pos.Z()*pos.Z());
    m_totalEnergy += energy;
    m_moment1Seed += energy * pos;
    m_moment2Seed += energy * pos2;
}


void McIntegratingHit::setEnergyItems(const Double_t& totE, const Double_t *energyArr,
                                      const TVector3& moment1, const TVector3& moment2) {

    // Purpose and Method:  Set all energy member variables

    m_totalEnergy = totE;
    m_energyArray[McIntegratingHit::PRIMARY] = energyArr[McIntegratingHit::PRIMARY];
    m_energyArray[McIntegratingHit::ELECTRON] = energyArr[McIntegratingHit::ELECTRON];
    m_energyArray[McIntegratingHit::POSITRON] = energyArr[McIntegratingHit::POSITRON];
    m_moment1Seed = moment1 * totE;
    m_moment2Seed = moment2 * totE;
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


const McXtalEnergyDep* McIntegratingHit::getXtalEnergyDep(const TString& key) const
{
    TObject* xtalEnergyDep = m_xtalEnergyMap.FindObject(key);
    return dynamic_cast<McXtalEnergyDep*>(xtalEnergyDep);
}

McXtalEnergyDep* McIntegratingHit::getXtalEnergyDep(const TString& key)
{
    TObject* xtalEnergyDep = m_xtalEnergyMap.FindObject(key);

    // If no XtalEnergyDep object, then make one
    if (!xtalEnergyDep)
    {
        TObjString* newKey = new TObjString(key);

        xtalEnergyDep = new McXtalEnergyDep();
        m_xtalEnergyMap.Add(newKey, xtalEnergyDep);
    }

    return dynamic_cast<McXtalEnergyDep*>(xtalEnergyDep);
}
    
/// Code for the XtalEnergyDep subclass
/// Retrieve the energy-weighted first moments of the position
const TVector3 McXtalEnergyDep::getMoment1() const
{
    return m_moment1seed * (1./m_totalEnergy);
}

TVector3 McXtalEnergyDep::getMoment1()
{
    return m_moment1seed * (1./m_totalEnergy);
}
    
/// Retrieve the energy-weighted second moments of the position
const TVector3 McXtalEnergyDep::getMoment2() const
{
    return m_moment2seed * (1./m_totalEnergy);
}

TVector3 McXtalEnergyDep::getMoment2()
{
    return m_moment2seed * (1./m_totalEnergy);
}

// add energy/position for this diode
void McXtalEnergyDep::addEnergyItem(Double_t totalEnergy, Double_t directEnergy, const TVector3& position)
{
    /// Increment the energy accumulators
    m_totalEnergy  += totalEnergy;
    m_directEnergy += directEnergy;

    /// Get position squared for moments accumulators
    TVector3 position2 = TVector3(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());

    /// Accumulate the moments
    m_moment1seed += totalEnergy * position;
    m_moment2seed += totalEnergy * position2;
}
    
void McXtalEnergyDep::initialize(Double_t totalEnergy, Double_t directEnergy, const TVector3& mom1, const TVector3& mom2)
{
    m_totalEnergy  = totalEnergy;
    m_directEnergy = directEnergy;
    m_moment1seed  = mom1;
    m_moment2seed  = mom2;
}    

void McXtalEnergyDep::Clear( Option_t * )
{
    m_totalEnergy  = 0.;
    m_directEnergy = 0.;
    m_moment1seed  = TVector3(0.,0.,0.);
    m_moment2seed  = TVector3(0.,0.,0.);
}

// dummy data, just for tests
void McXtalEnergyDep::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

    Clear() ;
    
    double totE    = 5.5;
    double directE = 1.1;
    TVector3 moment1(1.0, 2.0, 3.0);
    TVector3 moment2(2.0, 4.0, 6.0);
    initialize(totE,directE,moment1,moment2) ;

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att,ref.get ## att,#att)

Bool_t McXtalEnergyDep::CompareInRange( const McXtalEnergyDep& ref, const std::string & name ) const {

    Bool_t result = true ;
    
    result = COMPARE_IN_RANGE(TotalEnergy()) && result ;
    
    result = COMPARE_IN_RANGE(DirectEnergy()) && result ;
    result = COMPARE_IN_RANGE(Moment1()) && result ;
    result = COMPARE_IN_RANGE(Moment2()) && result ;
      
    if (!result) 
    {
        if ( name == "" ) 
        {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else 
        {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }

    return result ;
}

void McXtalEnergyDep::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    UInt_t p = 2;
    cout.precision(p);
    cout << "Total Energy:  " << m_totalEnergy << endl;
    cout << "Direct Energy: " << m_directEnergy << endl;
    cout << "Mom1: (" << m_moment1seed.X() << "," << m_moment1seed.Y() << ","
        << m_moment1seed.Z() << ")   ";
    cout << "Mom2: (" << m_moment2seed.X() << "," << m_moment2seed.Y() << ","
        << m_moment2seed.Z() << ")" << endl;

} 
