#include "mcRootData/McPositionHit.h"
#include <iostream>

ClassImp(McPositionHit)

McPositionHit::McPositionHit() :
m_statusFlags(0), m_depositedEnergy(0.), 
m_particleEnergy(0.), m_timeOfFlight(0.)
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
    m_entry = TVector3(0., 0., 0.);
    m_exit = TVector3(0., 0., 0.);
    m_volumeId.Clear();
}

void McPositionHit::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    m_volumeId.Print(option);
    cout.precision(2);
    cout << "Flags: " << m_statusFlags 
        << "    Dep Energy: " << m_depositedEnergy
        << "    Part Energy: " << m_particleEnergy
        << "    TOF:  " << m_timeOfFlight << endl;
    cout << "Entry: (" << m_entry.X() << ","
        << m_entry.Y() << "," << m_entry.Z() << ")"
        << "   Exit: (" << m_exit.X() << "," 
        << m_exit.Y() << "," << m_exit.Z() << ")" << endl;
    cout << "McParticle Ref" << endl;
    if (getMcParticle()!= 0) getMcParticle()->Print();
    cout << "Origin Particle Ref" << endl;
    if (getOriginMcParticle() != 0) getOriginMcParticle()->Print();

}

Double_t McPositionHit::getDirectionCosine() const
{
    Double_t dx = m_exit.x()-m_entry.x();
    Double_t dy = m_exit.y()-m_entry.y();
    Double_t dz = m_exit.z()-m_entry.z();
    return (dz / sqrt(dx * dx + dy * dy));
}


void McPositionHit::initialize(Double_t edep, const VolumeIdentifier &id, 
                         const TVector3& entry, const TVector3& exit,
                         McParticle *mc, 
                         McParticle *origin, Double_t pE, Double_t tof, UInt_t flags)
{
    m_depositedEnergy = edep;
    m_volumeId = id;
    m_entry.SetXYZ(entry.X(), entry.Y(), entry.Z());
    m_exit.SetXYZ(exit.X(), exit.Y(), exit.Z());
    m_mcParticle = mc;
    m_originMcParticle = origin;
    m_particleEnergy = pE;
    m_timeOfFlight = tof;
    m_statusFlags = flags;
}

Bool_t McPositionHit::needDigi() const
{
    return m_statusFlags & McParticle::NEED_DIGI;
}
