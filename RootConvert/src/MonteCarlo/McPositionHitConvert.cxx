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

}
