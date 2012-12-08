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
                       const CLHEP::HepLorentzVector& initialMomentum,
                       const CLHEP::HepLorentzVector& finalMomentum,
                       const HepPoint3D& initialPosition,
                       const HepPoint3D& finalPosition, 
                       const std::string process)
{
    initialize(mother, id, flags, initialMomentum, initialPosition,process);
    finalize(finalMomentum, finalPosition);
}

void McParticle::initialize( McParticle* mother,         
                      StdHepId id, 
        unsigned int flags,
        const CLHEP::HepLorentzVector& initialMomentum,
        const HepPoint3D& initialPosition, const std::string process)
{
    m_mother = mother;
    m_particleID = id;
    m_statusFlags = flags;
    m_initialFourMomentum = initialMomentum;
    m_initialPosition = initialPosition;
    m_process = process;
    if( mother != this) mother->m_daughters.push_back(this);
}

void McParticle::finalize(const CLHEP::HepLorentzVector& finalMomentum,
        const HepPoint3D& finalPosition)
{
    m_finalFourMomentum = finalMomentum;
    m_finalPosition = finalPosition;
}


void McParticle::transform(
        const CLHEP::HepLorentzVector& initialMomentum,
        const CLHEP::HepLorentzVector& finalMomentum,
        const HepPoint3D& initialPosition,
        const HepPoint3D& finalPosition)
{
    m_initialFourMomentum = initialMomentum;
    m_initialPosition     = initialPosition;
    m_finalFourMomentum   = finalMomentum;
    m_finalPosition       = finalPosition;
}


void McParticle::setMother(const SmartRef<McParticle> m)
{
  m_mother = m;
}

const HepPoint3D& McParticle::initialPosition()const
{
    return m_initialPosition;
}
const HepPoint3D& McParticle::finalPosition()const
{
    return m_finalPosition;
}
const CLHEP::HepLorentzVector&  McParticle::initialFourMomentum()const
{
    return m_initialFourMomentum;
}
const CLHEP::HepLorentzVector&  McParticle::finalFourMomentum()const
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

/// Remove daughters when in prune mode
void McParticle::removeDaughter(const SmartRef<McParticle> mcPart)
{
    // Particle to delete was most likely last one added, use reverse iterator
    SmartRefVector<Event::McParticle>::reverse_iterator daughtIter;
    //std::vector<SmartRef<Event::McParticle> >::reverse_iterator daughtIter;
    for(daughtIter = m_daughters.rbegin();daughtIter != m_daughters.rend();daughtIter++)
    {
        if (mcPart == *daughtIter)
        {
            SmartRefVector<Event::McParticle>::iterator forwardIter = (++daughtIter).base();
            m_daughters.erase(forwardIter);
            break;
        }
    }
    return;
}


/// access to the list of daughters: null if none
const SmartRefVector<McParticle>& McParticle::daughterList()const
{
    return m_daughters;
}


} // namespace



