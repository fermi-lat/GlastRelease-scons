#include "mcRootData/McPositionHit.h"

ClassImp(McPositionHit)

McPositionHit::McPositionHit()
: m_depositedEnergy(0.), m_particleEnergy(0.), m_timeOfFlight(0.), m_statusFlags(0)
{
}

McPositionHit::~McPositionHit()
{
    Clear();
}

void McPositionHit::Clear(Option_t *option)
{
    m_depositedEnergy = 0.;
    m_particleEnergy = 0;
    m_timeOfFlight = 0.;
    m_statusFlags = 0;
}


Double_t McPositionHit::getDirectionCosine() const
{
    Double_t dx = m_exit.x()-m_entry.x();
    Double_t dy = m_exit.y()-m_entry.y();
    Double_t dz = m_exit.z()-m_entry.z();
    return (dz / sqrt(dx * dx + dy * dy));
}


void McPositionHit::initialize(Double_t edep, VolumeIdentifier id, 
                         const TVector3& entry, const TVector3& exit,
                         McParticle *mc, 
                         McParticle *origin, Double_t pE, Double_t tof, UInt_t flags)
{
    m_depositedEnergy = edep;
    m_volumeId = id;
    m_entry = entry;
    m_exit = exit;
    m_mcParticle = mc;
    m_originMcParticle = origin;
    m_particleEnergy = pE;
    m_timeOfFlight = tof;
    m_statusFlags = flags;
}

Bool_t McPositionHit::needDigi() const
{
    return m_statusFlags & McParticle::StatusBits::NEED_DIGI;
}
