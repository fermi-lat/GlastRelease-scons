//                                                                       
// The McParticle class provides an interface to a particle for use in
// GLAST Monte Carlo applications
//
// $Header$

#include "mcRootData/McParticle.h"

ClassImp(McParticle)

McParticle::McParticle() : 
m_particleId(0), m_statusFlags(0)
{
    Clear();
}

McParticle::McParticle(const McParticle &p) {
    // Purpose and Method:  The copy constructor.  This will make it easier to
    //    read in the McParticles within the simulation.

    m_particleId = p.m_particleId;
    m_statusFlags = p.m_statusFlags;
    m_finalPosition = TVector3(p.m_finalPosition.X(), p.m_finalPosition.Y(),
        p.m_finalPosition.Z());
    m_initialFourMomentum = TLorentzVector(p.m_initialFourMomentum.X(),
        p.m_initialFourMomentum.Y(), p.m_initialFourMomentum.Z(), 
        p.m_initialFourMomentum.T());
    m_finalFourMomentum = TLorentzVector(p.m_finalFourMomentum.X(),
        p.m_finalFourMomentum.Y(), p.m_finalFourMomentum.Z(),
        p.m_finalFourMomentum.T());
    m_mother = p.m_mother;
    TRefArrayIter daughterIter(&p.m_daughters);
    McParticle *curRef;
    while (curRef = (McParticle*)daughterIter.Next()){
        m_daughters.Add(curRef);
    }
}

McParticle::~McParticle() {
    Clear();
}

void McParticle::Clear(Option_t *option) {
    m_particleId = 0;
    m_statusFlags = 0;
    m_finalPosition = TVector3(0., 0., 0.);
    m_initialFourMomentum = TLorentzVector(0., 0., 0., 0.);
    m_finalFourMomentum = TLorentzVector(0., 0., 0., 0.);
    m_daughters.Clear();
}

void McParticle::initialize( McParticle* mother, Int_t id, UInt_t statusBits,
                      const TLorentzVector& initMom,
                      const TLorentzVector& finalMom,
                      const TVector3& finalPos)
{
    m_mother = mother;
    m_particleId = id;
    m_statusFlags = statusBits;
    m_initialFourMomentum.SetXYZT(initMom.X(), initMom.Y(), initMom.Z(), initMom.T());
    m_finalFourMomentum.SetXYZT(finalMom.X(), finalMom.Y(), finalMom.Z(), finalMom.T());
    m_finalPosition.SetXYZ(finalPos.X(), finalPos.Y(), finalPos.Z());
    if ( mother == 0 ) return;
    if( mother != this) mother->m_daughters.Add(this);
}


McParticle* McParticle::getMother() {

    return ( (McParticle*)(m_mother.GetObject()) );
}

const McParticle* McParticle::getDaughter(Int_t index) const {
    return ( (McParticle*)m_daughters.At(index) );
}

Int_t McParticle::getParticleProperty() const {
    return m_statusFlags;
}

Bool_t McParticle::primaryParticle() const
{
  return (m_statusFlags & PRIMARY)==PRIMARY;
}


const TVector3& McParticle::getInitialPosition()const
{
    return ((McParticle*)m_mother.GetObject())->m_finalPosition;
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
