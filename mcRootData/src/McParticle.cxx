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
}

McParticle::~McParticle() {
}

void McParticle::Clear(Option_t *option) {
    m_statusFlags = 0;
}

void McParticle::initialize( McParticle* mother, Int_t id, UInt_t statusBits,
                      const TLorentzVector& initalMomentum,
                      const TLorentzVector& finalMomentum,
                      const TVector3& finalPosition)
{
    m_mother = mother;
    m_particleId = id;
    m_statusFlags = statusBits;
    m_initialFourMomentum = initalMomentum;
    m_finalFourMomentum = finalMomentum;
    m_finalPosition = finalPosition;
    if( mother != this) mother->m_daughters.Add(this);
}


const McParticle* McParticle::getMother() const {

    return ( (McParticle*)m_mother.GetObject() );
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

