// $Header$

#include <iostream>
#include "Event/MonteCarlo/McParticle.h"
#include "Event/Utilities/CLHEPStreams.h"


namespace Event {

/// Retrieve particle property
McParticle::StdHepId McParticle::particleProperty() const
{
  return m_particleID;
}


/// Retrieve whether this is a primary particle
bool McParticle::primaryParticle() const
{
  return (m_statusFlags & PRIMARY)==PRIMARY;
}

void McParticle::init( McParticle* mother,         
                      StdHepId id, 
        unsigned int flags,
        const HepLorentzVector& initialMomentum,
        const HepLorentzVector& finalMomentum,
        const HepPoint3D& initialPosition,
        const HepPoint3D& finalPosition)
{
    initialize(mother, id, flags, initialMomentum, initialPosition);
    finalize(finalMomentum, finalPosition);
}

void McParticle::initialize( McParticle* mother,         
                      StdHepId id, 
        unsigned int flags,
        const HepLorentzVector& initialMomentum,
        const HepPoint3D& initialPosition)
{
    m_mother = mother;
    m_particleID = id;
    m_statusFlags = flags;
    m_initialFourMomentum = initialMomentum;
    m_initialPosition = initialPosition;
    if( mother != this) mother->m_daughters.push_back(this);
}

void McParticle::finalize(const HepLorentzVector& finalMomentum,
        const HepPoint3D& finalPosition)
{
    m_finalFourMomentum = finalMomentum;
    m_finalPosition = finalPosition;
}


const HepPoint3D& McParticle::initialPosition()const
{
    return m_initialPosition;
}
const HepPoint3D& McParticle::finalPosition()const
{
    return m_finalPosition;
}
const HepLorentzVector&  McParticle::initialFourMomentum()const
{
    return m_initialFourMomentum;
}
const HepLorentzVector&  McParticle::finalFourMomentum()const
{
    return m_finalFourMomentum;
}

unsigned int McParticle::statusFlags()const{
    return m_statusFlags;
}


/// access to the mother particle
const McParticle& McParticle::mother()const
{
    return *m_mother;
}


/// access to the list of daughters: null if none
const SmartRefVector<McParticle>& McParticle::daughterList()const
{
    return m_daughters;
}


} // namespace



