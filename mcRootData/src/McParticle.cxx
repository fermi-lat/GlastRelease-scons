//                                                                       
// The McParticle class provides an interface to a particle for use in
// GLAST Monte Carlo applications
//
// $Header$

#include "mcRootData/McParticle.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

#include "McObjectManager.h"

ClassImp(McParticle)

McParticle::McParticle() : 
m_particleId(0), m_statusFlags(0)
{
    Clear();
}

McParticle::McParticle(const McParticle &p)
 : TObject(p) {
    // Purpose and Method:  The copy constructor.  This will make it easier to
    //    read in the McParticles within the simulation.

    m_particleId = p.m_particleId;
    m_statusFlags = p.m_statusFlags;
    m_initialPosition = p.m_initialPosition;
    m_finalPosition = p.m_finalPosition;
    m_initialFourMomentum = p.m_initialFourMomentum;
    m_finalFourMomentum = p.m_finalFourMomentum;
    m_mother = p.m_mother;
    m_process = p.m_process;
    // TRefArrayIter::Next does not seem to be implemented in Root 3.02.03
    // or ROOT 3.03.04 - just iterating over the entries for now.
    //TRefArrayIter daughterIter(&p.m_daughters);
    McParticle *curRef;
    Int_t iPart;
    for (iPart = 0; iPart < p.m_daughters.GetEntries(); iPart++) {
    //while ( (curRef = (McParticle*)daughterIter.Next()) ){
        curRef = (McParticle*)p.m_daughters.At(iPart);
        m_daughters.Add(curRef);
    }
}

McParticle::~McParticle() {
    Clear();
}

void* McParticle::operator new(size_t size)
{
    McParticle* temp = McObjectManager::getPointer()->getNewMcParticle();

    // Since we recycle, make sure these member functions are cleared
    temp->m_daughters.Clear();
    temp->m_process = "";

    return temp;
}

void* McParticle::operator new(size_t size, void* vp)
{
    return vp;
}

McParticle& McParticle::operator =(const McParticle& rhs)
{
    m_particleId           = rhs.m_particleId;
    m_statusFlags          = rhs.m_statusFlags;
    m_initialPosition      = rhs.m_initialPosition;
    m_finalPosition        = rhs.m_finalPosition;
    m_initialFourMomentum  = rhs.m_initialFourMomentum;
    m_finalFourMomentum    = rhs.m_finalFourMomentum;
    m_mother               = rhs.m_mother;
    m_process              = rhs.m_process; 

    for(int idx = 0; idx < rhs.m_daughters.GetEntries(); idx++)
        m_daughters.Add(rhs.m_daughters[idx]);

    return *this;
}

void McParticle::Clear(Option_t *) {
    m_particleId = 0;
    m_statusFlags = 0;
    m_initialPosition = TVector3(0., 0., 0.);
    m_finalPosition = TVector3(0., 0., 0.);
    m_initialFourMomentum = TLorentzVector(0., 0., 0., 0.);
    m_finalFourMomentum = TLorentzVector(0., 0., 0., 0.);
    m_daughters.Clear();
}

// dummy data, just for tests
void McParticle::Fake( Int_t ievent, UInt_t rank, Float_t randNum ) {

    Clear() ;
    Float_t f = Float_t(rank);
    Float_t fr = f*randNum;
    TLorentzVector initMom(f, f, f, f);
    TLorentzVector finalMom(fr+ievent, fr+ievent, fr+ievent, fr+ievent);
    TVector3 initPos(randNum, f*2.0*randNum, f*4.0*randNum);
    TVector3 finalPos(f*randNum, f*randNum, f*randNum);
    initialize(this, rank, 0, initMom, finalMom, initPos, finalPos);

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t McParticle::CompareInRange( const McParticle & ref, const std::string & name ) const {

    Bool_t result = true ;
    
    result = COMPARE_IN_RANGE(ParticleId) && result ;
    result = COMPARE_IN_RANGE(StatusFlags) && result ;
    
    result = COMPARE_IN_RANGE(InitialFourMomentum) && result ;
    result = COMPARE_IN_RANGE(FinalFourMomentum) && result ;

    result = COMPARE_IN_RANGE(InitialPosition) && result ;
    result = COMPARE_IN_RANGE(FinalPosition) && result ;
        
    result = COMPARE_IN_RANGE(Process) && result ;
    
    // for fake data and comparison
    // particles should point to themselves
    rootdatautil::CompareInRange(ref.GetUniqueID(),ref.getMother()->GetUniqueID(),"Ref Mother") ;
    rootdatautil::CompareInRange(GetUniqueID(),getMother()->GetUniqueID(),"Mother") ;
      
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

void McParticle::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "ParticleId: " << m_particleId;
    cout << " StatusFlag: " << m_statusFlags << endl;
    cout << "Process: " << m_process.Data() << endl;
    cout << "InitialPos: (" << m_initialPosition.X() << "," << m_initialPosition.Y()
        << "," << m_initialPosition.Z() << ") ";
    cout << "FinalPos: (" << m_finalPosition.X() << "," << m_finalPosition.Y()
        << "," << m_finalPosition.Z() << ")" << endl;
    cout << "InitMom: (" << m_initialFourMomentum.X() << "," << m_initialFourMomentum.Y() << ","
        << m_initialFourMomentum.Z() << "," << m_initialFourMomentum.T() << ")    ";
    cout << "FinalMom: (" << m_finalFourMomentum.X() << "," << m_finalFourMomentum.Y() << ","
        << m_finalFourMomentum.Z() << "," << m_finalFourMomentum.T() << ")" << endl;
    cout << "Mother:  ";
    m_mother.Print();
    cout << "Number of Daughters: " << m_daughters.GetEntries() << endl;
}


void McParticle::initialize( McParticle* mother, Int_t id, UInt_t statusBits,
                      const TLorentzVector& initMom,
                      const TLorentzVector& finalMom,
                      const TVector3& initPos,
                      const TVector3& finalPos,
                      const char* process)
{
    m_mother = mother;
    m_particleId = id;
    m_statusFlags = statusBits;
    m_initialFourMomentum = initMom;
    m_finalFourMomentum = finalMom;
    m_initialPosition = initPos;
    m_finalPosition = finalPos;
    m_process = TString(process);
    if ( mother == 0 ) return;
    if( mother != this) mother->m_daughters.Add(this);
}


const McParticle* McParticle::getMother() const {

    return ( (McParticle*)(m_mother.GetObject()) );
}

const McParticle* McParticle::getDaughter(Int_t index) const {
    return ( (McParticle*)m_daughters.At(index) );
}


Bool_t McParticle::primaryParticle() const
{
  return (m_statusFlags & PRIMARY)==PRIMARY;
}


const TVector3& McParticle::getInitialPosition()const
{
    return m_initialPosition;
}


const TVector3& McParticle::getFinalPosition() const {
    return m_finalPosition;
}

    
const TLorentzVector& McParticle::getInitialFourMomentum() const { 
    return m_initialFourMomentum; 
}

const TLorentzVector& McParticle::getFinalFourMomentum() const {
    return m_finalFourMomentum;
};

const TString& McParticle::getProcess() const {
    return m_process;
}

