// $Header$

#include <iostream>
#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"


namespace mc {

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
        const HepPoint3D& finalPosition)
{
    m_mother = mother;
    m_particleID = id;
    m_statusFlags = flags;
    m_initialFourMomentum = initialMomentum;
    m_finalFourMomentum = finalMomentum;
    m_finalPosition = finalPosition;
    if( mother != this) mother->m_daughters.push_back(this);
}

const HepPoint3D& McParticle::initialPosition()const
{
    return m_mother->m_finalPosition;
}
const HepPoint3D& McParticle::finalPosition()const
{
    return m_finalPosition;
}
const HepLorentzVector&  McParticle::initialFourMomemtum()const
{
    return m_initialFourMomentum;
}
const HepLorentzVector&  McParticle::finalFourMomemtum()const
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
const SmartRefVector<McParticle>& McParticle::daugherList()const
{
    return m_daughters;
}


} // namespace



