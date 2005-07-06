#include <RootConvert/MonteCarlo/McPositionHitConvert.h>
#include <RootConvert/Utilities/Toolkit.h>
#include "TRandom.h"

namespace RootPersistence {

void convert( const Event::McPositionHit & tdsHit, McPositionHit & rootHit ) {

    //bool primaryOrigin = tdsPositionHit.primaryOrigin();
    //log << MSG::INFO << "primaryOrigin " << primaryOrigin << endreq;
        
    //bool caloShowerOrigin = tdsPositionHit.caloShowerOrigin();
    //log << MSG::INFO << "calShowerO " << caloShowerOrigin << endreq;
        
    //bool needDigi = tdsPositionHit.needDigi();
    //log << MSG::INFO << "need Digi " << needDigi << endreq;
               
    // not to be managed by RootConvert
    McParticle * mcPartRoot = 0 ;
    McParticle * originRoot = 0 ;
    UInt_t flagsRoot = 0 ;

    rootHit.initialize(
        tdsHit.getMcParticleId(),
        tdsHit.getOriginMcParticleId(),
        tdsHit.depositedEnergy(),
        convert(tdsHit.volumeID()), 
        convert(tdsHit.entryPoint()),
        convert(tdsHit.exitPoint()),
        convert(tdsHit.globalEntryPoint()),
        convert(tdsHit.globalExitPoint()), 
        mcPartRoot,originRoot,
        convert(tdsHit.particleFourMomentum()),
        tdsHit.timeOfFlight(),
        flagsRoot) ;

}
 
void convert( const McPositionHit & rootHit, Event::McPositionHit & tdsHit ) {

    tdsHit.init(
        rootHit.getDepositedEnergy(),
        convert(rootHit.getVolumeId()),
        convert(rootHit.getEntryPosition()),
        convert(rootHit.getExitPosition()),
        convert(rootHit.getGlobalEntryPosition()),
        convert(rootHit.getGlobalExitPosition()) ) ;

    tdsHit.setParticle4Momentum(convert(rootHit.getParticle4Momentum()));
    tdsHit.setTimeOfFlight(rootHit.getTimeOfFlight());
    tdsHit.setMcParticleId(rootHit.getMcParticleId()) ;
    tdsHit.setOriginMcParticleId(rootHit.getOriginMcParticleId()) ;

}

void random( McPositionHit & hit, UInt_t ipart, Float_t randNum ) {

    Float_t f = Float_t(ipart);
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
    
    hit.initialize(particleId,originId,depE,id, 
      entry,exit,gEntry,gExit,
      mcPart,originMcPart,part4Mom,tof,flags) ;

}

bool compare( const McPositionHit & hit1, const McPositionHit & hit2 ) {

    bool result = true ;
    
    result = result && compare(hit1.getMcParticleId(),hit2.getMcParticleId(),"McParticleId") ;
    result = result && compare(hit1.getOriginMcParticleId(),hit2.getOriginMcParticleId(),"OriginMcParticleId") ;

    VolumeIdentifier id1 = hit1.getVolumeId() ;
    VolumeIdentifier id2 = hit2.getVolumeId() ;
    result = result && compare(id1.name(),id2.name(),"VolumeId Name") ;
    result = result && compare(id1.size(),id2.size(),"VolumeId Size") ;
    result = result && compare(id1.getBits0to31(),id2.getBits0to31(),"VolumeId Bits0to31") ;
    result = result && compare(id1.getBits32to63(),id2.getBits32to63(),"VolumeId Bits32to63") ;
    
    result = result && compare(hit1.getEntryPosition(),hit2.getEntryPosition(),"Entry") ;
    result = result && compare(hit1.getExitPosition(),hit2.getExitPosition(),"Exit") ;
    result = result && compare(hit1.getGlobalEntryPosition(),hit2.getGlobalEntryPosition(),"GlobalEntry") ;
    result = result && compare(hit1.getGlobalExitPosition(),hit2.getGlobalExitPosition(),"GlobalExit") ;

    result = result && compare(hit1.getDepositedEnergy(),hit2.getDepositedEnergy(),"DepositedEnergy") ;
    result = result && compare(hit1.getParticleEnergy(),hit2.getParticleEnergy(),"ParticleEnergy") ;
    result = result && compare(hit1.getTimeOfFlight(),hit2.getTimeOfFlight(),"TimeOfFlight") ;

    return result ;
}

}
