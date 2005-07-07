#include <mcRootData/McPositionHit.h>
#include <commonRootData/RootDataUtil.h>
#include <iostream>

ClassImp(McPositionHit)

McPositionHit::McPositionHit() :
m_depositedEnergy(0.), m_timeOfFlight(0.), 
m_particleFourMomentum(0.,0.,0.,0.), m_statusFlags(0)
{
}

McPositionHit::~McPositionHit()
{
    Clear();
}

void McPositionHit::Clear(Option_t *)
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

// dummy data, just for tests
void McPositionHit::Fake( UInt_t rank, Float_t randNum ) {

    Float_t f = Float_t(rank);
    Float_t fr = f*randNum;
    TVector3 entry(1., 1., 1.);
    TVector3 exit(fr, fr, fr);
    TVector3 gEntry(3., 3., 3.);
    TVector3 gExit(fr*2., fr*2., fr*2.);
            
    VolumeIdentifier id ;
    id.append(1) ;
            
    Double_t depE = randNum;
    Double_t partE = randNum*0.1 ;
    TVector3 partMom = entry-exit ;
    TLorentzVector part4Mom(partMom,partE) ;
    Double_t tof = randNum*0.4 ;
    UInt_t flags = 0 ;
    McParticle * mcPart = 0 ;
    McParticle * originMcPart = 0 ;
    Int_t particleId = 7 ;
    Int_t originId = -13 ;
    
    initialize(particleId,originId,depE,id, 
      entry,exit,gEntry,gExit,
      mcPart,originMcPart,part4Mom,tof,flags) ;

}

Bool_t McPositionHit::CompareInRange( const McPositionHit & hit ) const {

    Bool_t result = true ;
    
    result = result && rootdatautil::CompareInRange(getMcParticleId(),hit.getMcParticleId(),"McParticleId") ;
    result = result && rootdatautil::CompareInRange(getOriginMcParticleId(),hit.getOriginMcParticleId(),"OriginMcParticleId") ;

    VolumeIdentifier id1 = getVolumeId() ;
    VolumeIdentifier id2 = hit.getVolumeId() ;
    result = result && rootdatautil::CompareInRange(id1.name(),id2.name(),"VolumeId Name") ;
    result = result && rootdatautil::CompareInRange(id1.size(),id2.size(),"VolumeId Size") ;
    result = result && rootdatautil::CompareInRange(id1.getBits0to31(),id2.getBits0to31(),"VolumeId Bits0to31") ;
    result = result && rootdatautil::CompareInRange(id1.getBits32to63(),id2.getBits32to63(),"VolumeId Bits32to63") ;
    
    result = result && rootdatautil::CompareInRange(getEntryPosition(),hit.getEntryPosition(),"Entry") ;
    result = result && rootdatautil::CompareInRange(getExitPosition(),hit.getExitPosition(),"Exit") ;
    result = result && rootdatautil::CompareInRange(getGlobalEntryPosition(),hit.getGlobalEntryPosition(),"GlobalEntry") ;
    result = result && rootdatautil::CompareInRange(getGlobalExitPosition(),hit.getGlobalExitPosition(),"GlobalExit") ;

    result = result && rootdatautil::CompareInRange(getDepositedEnergy(),hit.getDepositedEnergy(),"DepositedEnergy") ;
    result = result && rootdatautil::CompareInRange(getParticleEnergy(),hit.getParticleEnergy(),"ParticleEnergy") ;
    result = result && rootdatautil::CompareInRange(getTimeOfFlight(),hit.getTimeOfFlight(),"TimeOfFlight") ;

    if (!result) {
        std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
    }
    return result ;

}

void McPositionHit::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    m_volumeId.Print(option);
    std::cout.precision(2);
    std::cout << "Flags: " << m_statusFlags 
        << "    Dep Energy: " << m_depositedEnergy
        << "    Part Energy: " << m_particleFourMomentum.E()
        << "    TOF:  " << m_timeOfFlight << std::endl;
    std::cout << "Local Entry: (" << m_entry.X() << ","
        << m_entry.Y() << "," << m_entry.Z() << ")"
        << "  Local Exit: (" << m_exit.X() << "," 
        << m_exit.Y() << "," << m_exit.Z() << ") Direction Cosine: " 
        << getDirectionCosine() << std::endl;
    std::cout << "Global Entry: (" << m_globalEntry.X() << ","
        << m_globalEntry.Y() << "," << m_globalEntry.Z() << ")"
        << "  Global Entry: (" << m_globalExit.X() << ","
        << m_globalExit.Y() << "," << m_globalExit.Z() << ")" << std::endl;
    std::cout << "McParticleId " << m_mcParticleId << std::endl;
    if (getMcParticle() != 0) {
        cout << "McParticle Ref" << std::endl;
        getMcParticle()->Print();
    }
    std::cout << "Origin ParticleId " << m_originMcParticleId << endl;
    if (getOriginMcParticle() != 0) {
        cout << "Origin Particle Ref" << std::endl;
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
                               const TLorentzVector& p4Mom, 
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
                               McParticle *mc, McParticle *origin, const TLorentzVector& p4Mom,
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
