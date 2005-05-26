#include "mcRootData/McPositionHit.h"
#include <iostream>

ClassImp(McPositionHit)

McPositionHit::McPositionHit() :
m_depositedEnergy(0.), m_particleFourMomentum(0.,0.,0.,0.), 
m_timeOfFlight(0.), m_statusFlags(0)
{
}

McPositionHit::~McPositionHit()
{
    Clear();
}

void McPositionHit::Clear(Option_t *option)
{
    m_depositedEnergy = 0.;
    m_particleFourMomentum = TLorentzVector(0.,0.,0.,0.);
    m_timeOfFlight = 0.;
    m_statusFlags = 0;
    m_entry = TVector3(0., 0., 0.);
    m_exit = TVector3(0., 0., 0.);
    m_globalEntry = TVector3(0., 0., 0.);
    m_globalExit = TVector3(0., 0., 0.);
    m_volumeId.Clear();
}

void McPositionHit::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    m_volumeId.Print(option);
    cout.precision(2);
    cout << "Flags: " << m_statusFlags 
        << "    Dep Energy: " << m_depositedEnergy
        << "    Part Energy: " << m_particleFourMomentum.E()
        << "    TOF:  " << m_timeOfFlight << endl;
    cout << "Local Entry: (" << m_entry.X() << ","
        << m_entry.Y() << "," << m_entry.Z() << ")"
        << "  Local Exit: (" << m_exit.X() << "," 
        << m_exit.Y() << "," << m_exit.Z() << ") Direction Cosine: " 
        << getDirectionCosine() << endl;
    cout << "Global Entry: (" << m_globalEntry.X() << ","
        << m_globalEntry.Y() << "," << m_globalEntry.Z() << ")"
        << "  Global Entry: (" << m_globalExit.X() << ","
        << m_globalExit.Y() << "," << m_globalExit.Z() << ")" << endl;
    cout << "McParticleId " << m_mcParticleId << endl;
    if (getMcParticle() != 0) {
        cout << "McParticle Ref" << endl;
        getMcParticle()->Print();
    }
    cout << "Origin ParticleId " << m_originMcParticleId << endl;
    if (getOriginMcParticle() != 0) {
        cout << "Origin Particle Ref" << endl;
        getOriginMcParticle()->Print();
    }

}

Double_t McPositionHit::getDirectionCosine() const
{ 
    TVector3 dir = (m_exit - m_entry).Unit();
    return dir.Z();
}


void McPositionHit::initialize(Int_t particleId, Double_t edep, 
                               const VolumeIdentifier &volId,
                               const TVector3& entry, const TVector3& exit,
                               McParticle *mc, McParticle *origin, 
                               TLorentzVector& p4Mom, 
                               Double_t tof, UInt_t flags)
{
    m_mcParticleId = particleId;
    m_depositedEnergy = edep;
    m_volumeId = volId;
    m_entry = entry;
    m_exit = exit;
    m_mcParticle = mc;
    m_originMcParticle = origin;
    m_particleFourMomentum = p4Mom;
    m_timeOfFlight = tof;
    m_statusFlags = flags;
}

void McPositionHit::initialize(Int_t mcParticleId, Int_t originParticleId, 
                               Double_t edep, const VolumeIdentifier& volId, 
                               const TVector3& entry, const TVector3& exit,
                               const TVector3& gEntry, const TVector3& gExit, 
                               McParticle *mc, McParticle *origin, TLorentzVector& p4Mom,
                               Double_t tof, UInt_t flags)
{
    m_mcParticleId = mcParticleId;
    m_originMcParticleId = originParticleId;
    m_depositedEnergy = edep;
    m_volumeId = volId;
    m_entry = entry;
    m_exit = exit;
    m_globalEntry = gEntry;
    m_globalExit = gExit;
    m_mcParticle = mc;
    m_originMcParticle = origin;
    m_particleFourMomentum = p4Mom;
    m_timeOfFlight = tof;
    m_statusFlags = flags;
}


Bool_t McPositionHit::needDigi() const
{
    return m_statusFlags & McParticle::NEED_DIGI;
}
