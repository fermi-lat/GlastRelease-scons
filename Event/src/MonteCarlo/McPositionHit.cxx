// $Header$

#include <iostream>
#include <math.h>
#include "GlastEvent/MonteCarlo/McPositionHit.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"


// FIXME!!:
// The next `using' directive is an ad-hoc declaration for the transition
// phase to the namespace `GlastEvent'.
// After the transition is completed, it should be removed and each
// function definition should have the namespace identifier.
using namespace GlastEvent;

/// Retrieve hit's direction cosine
double McPositionHit::directionCosine() const
{
  double dx = m_exit.x()-m_entry.x();
  double dy = m_exit.y()-m_entry.y();
  double dz = m_exit.z()-m_entry.z();
  return dz / sqrt(dx * dx + dy * dy);
}





const VolumeID McPositionHit::volumeID() const
{
    return m_volumeID;
}
/// Update cell identifier
void McPositionHit::setVolumeID( VolumeID value )
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
    return m_particleEnergy;
}

/// Set depositing particle's energy
void McPositionHit::setParticleEnergy(double value)
{
    m_particleEnergy = value;
}

/// Retrieve primary-origin flag
bool McPositionHit::primaryOrigin() const
{
    using GlastEvent::McConstants::ORIGIN_PRIMARY;
    return m_packedFlags & ORIGIN_PRIMARY;
}
/// Update primary-origin flag
void McPositionHit::setPrimaryOrigin( bool value )
{
    using GlastEvent::McConstants::ORIGIN_PRIMARY;
    if (value){
        m_packedFlags |= ORIGIN_PRIMARY;
    } else {
        m_packedFlags &= ~ORIGIN_PRIMARY;
    }
}
/// Retrieve calorimeter-shower-origin flag
bool McPositionHit::caloShowerOrigin() const
{
    using GlastEvent::McConstants::ORIGIN_CALOSHOWER;
    return m_packedFlags & ORIGIN_CALOSHOWER;
}
/// Update calorimeter-shower-origin flag
void McPositionHit::setCaloShowerOrigin( bool value )
{
    using GlastEvent::McConstants::ORIGIN_CALOSHOWER;
    if (value){
        m_packedFlags |= ORIGIN_CALOSHOWER;
    } else {
        m_packedFlags &= ~ORIGIN_CALOSHOWER;
    }
}

/// Retrieve whether this hit should be digitized
bool McPositionHit::needDigi() const
{
    using GlastEvent::McConstants::NEED_DIGI;
    return m_packedFlags & NEED_DIGI;
}
/// Update whether this hit should be digitized
void McPositionHit::setNeedDigi( bool value )
{
    using GlastEvent::McConstants::NEED_DIGI;
    if (value){
        m_packedFlags |= NEED_DIGI;
    } else {
        m_packedFlags &= ~NEED_DIGI;
    }
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
