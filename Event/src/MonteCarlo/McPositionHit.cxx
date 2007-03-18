// $Header$

#include <iostream>
#include <math.h>
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/Utilities/CLHEPStreams.h"

namespace Event{

/// Retrieve hit's direction cosine
double McPositionHit::directionCosine() const
{
  double dx = m_exit.x()-m_entry.x();
  double dy = m_exit.y()-m_entry.y();
  double dz = m_exit.z()-m_entry.z();
  return dz / sqrt(dx * dx + dy * dy);
}


idents::VolumeIdentifier McPositionHit::volumeID() const
{
    return m_volumeID;
}
/// Update cell identifier
void McPositionHit::setVolumeID( idents::VolumeIdentifier value )
{
    m_volumeID = value;
}


/// Retrieve entry member
const HepPoint3D& McPositionHit::entryPoint() const
{
    return m_entry;
}
HepPoint3D& McPositionHit::entryPoint()
{
    return m_entry;
}


/// Update Entry member
void McPositionHit::setEntryPoint( const HepPoint3D& value )
{
    m_entry = value;
}


/// Retrieve exit point
const HepPoint3D& McPositionHit::exitPoint() const
{
    return m_exit;
}
HepPoint3D& McPositionHit::exitPoint()
{
    return m_exit;
}


/// Update exit point
void McPositionHit::setExitPoint( const HepPoint3D& value )
{
    m_exit = value;
}


/// Retrieve entry member in global coordinates
const HepPoint3D& McPositionHit::globalEntryPoint() const
{
    return m_globalEntry;
}
HepPoint3D& McPositionHit::globalEntryPoint()
{
    return m_globalEntry;
}


/// Update Entry member in global coordinates
void McPositionHit::setGlobalEntryPoint( const HepPoint3D& value )
{
    m_globalEntry = value;
}


/// Retrieve exit point in global coordinates
const HepPoint3D& McPositionHit::globalExitPoint() const
{
    return m_globalExit;
}
HepPoint3D& McPositionHit::globalExitPoint()
{
    return m_globalExit;
}


/// Update exit point in global coordinates
void McPositionHit::setGlobalExitPoint( const HepPoint3D& value )
{
    m_globalExit = value;
}


/// Retrieve deposited energy
double McPositionHit::depositedEnergy() const
{
    return m_depositedEnergy;
}


/// Update deposited energy
void McPositionHit::setDepositedEnergy( double value )
{
    m_depositedEnergy = value;
}

/// Retrieve depositing particle's energy
double McPositionHit::particleEnergy() const
{
    return m_particleFourMomentum.e();
}

/// Retrieve depositing particle's momentum
CLHEP::Hep3Vector McPositionHit::particleMomentum() const
{
    return m_particleFourMomentum.vect();
}

/// Retrieve depositing particle's four momentum
CLHEP::HepLorentzVector McPositionHit::particleFourMomentum() const
{
    return m_particleFourMomentum;
}

/// Set depositing particle's energy
//void McPositionHit::setParticleEnergy(double value)
//{
//    m_particleEnergy = value;
//}
void McPositionHit::setParticle4Momentum( const CLHEP::HepLorentzVector& fourMom)
{
    m_particleFourMomentum = fourMom;
    m_particleEnergy       = fourMom.e();
}

/// Retrieve member TOF
double McPositionHit::timeOfFlight() const
{
    return m_timeOfFlight;
}


/// Update TOF member
void McPositionHit::setTimeOfFlight( double value )
{
    m_timeOfFlight = value;
}


/// Retrieve pointer to McParticle (const or non-const)
const McParticle* McPositionHit::mcParticle() const
{
    return m_mcParticle; 
}
McParticle* McPositionHit::mcParticle()
{
    return m_mcParticle; 
}
/// Update pointer to McParticle (by a C++ pointer or a smart reference)
void McPositionHit::setMcParticle( McParticle* value )
{
    m_mcParticle = value; 
}
void McPositionHit::setMcParticle( SmartRef<McParticle> value )
{ 
    m_mcParticle = value; 
}


/// Retrieve pointer to the ancestor McParticle (const or non-const)
const McParticle* McPositionHit::originMcParticle() const
{
    return m_originMcParticle;
}
McParticle* McPositionHit::originMcParticle()
{
    return m_originMcParticle;
}
/// Update pointer to McParticle (by a C++ pointer or a smart reference)
void McPositionHit::setOriginMcParticle( McParticle* value )
{
    m_originMcParticle = value;
}
void McPositionHit::setOriginMcParticle( SmartRef<McParticle> value )
{
    m_originMcParticle = value;
}

}
