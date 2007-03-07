// File and Version Information:
//      $Header$
//
// Description:
// The McEvent class contains general information about a Monte Carlo
// Event, (Run/Event numbers, etc) as well as a list of McParticle
// instances generated during the simulated event.
//                                                
// Authors:  Dan Flath and Heather Kelly                       

#include <mcRootData/McEvent.h>
#include <commonRootData/RootDataUtil.h>
#include "TString.h"
#include "Riostream.h"

#include "McObjectManager.h"

// Allocate the TObjArrays just once
TObjArray *McEvent::m_staticParticleCol = 0;
TObjArray *McEvent::m_staticPositionHitCol = 0;
TObjArray *McEvent::m_staticIntegratingHitCol = 0;
TObjArray *McEvent::m_staticTrajectoryCol = 0;
TObjArray *McEvent::m_staticTkrStripCol = 0;

static const TString sourceUnknownStr = "NotStored";

ClassImp(McEvent)

McEvent::McEvent() : 
m_eventId(0), m_runId(0), m_sourceName(sourceUnknownStr)
{
    if (!m_staticParticleCol) m_staticParticleCol = new TObjArray();
    m_particleCol = m_staticParticleCol;
    
    if (!m_staticPositionHitCol) m_staticPositionHitCol = new TObjArray();
    m_positionHitCol = m_staticPositionHitCol;
    
    if (!m_staticTkrStripCol) m_staticTkrStripCol = new TObjArray();
    m_tkrStripCol = m_staticTkrStripCol;
    
    if (!m_staticIntegratingHitCol) m_staticIntegratingHitCol = new TObjArray();
    m_integratingHitCol = m_staticIntegratingHitCol;
    
    if (!m_staticTrajectoryCol) m_staticTrajectoryCol = new TObjArray();
    m_trajectoryCol = m_staticTrajectoryCol;
}

McEvent::~McEvent() {
    
    if (m_particleCol == m_staticParticleCol) m_staticParticleCol = 0;
    m_particleCol->Clear();  // We don't own these objects
    delete m_particleCol;
    m_particleCol = 0;

    if (m_positionHitCol == m_staticPositionHitCol) m_staticPositionHitCol = 0;
    m_positionHitCol->Clear();
    delete m_positionHitCol;
    m_positionHitCol = 0;

    if (m_tkrStripCol == m_staticTkrStripCol) m_staticTkrStripCol = 0;
    m_tkrStripCol->Clear();
    delete m_tkrStripCol;
    m_tkrStripCol = 0;

    if(m_integratingHitCol == m_staticIntegratingHitCol) m_staticIntegratingHitCol = 0;
    m_integratingHitCol->Clear();
    delete m_integratingHitCol;
    m_integratingHitCol = 0;

    if (m_trajectoryCol == m_staticTrajectoryCol) m_staticTrajectoryCol = 0;
    m_trajectoryCol->Clear();
    delete m_trajectoryCol;
    m_trajectoryCol = 0;

    m_sourceName = sourceUnknownStr;

    delete McObjectManager::getPointer();
}

void McEvent::initialize(UInt_t nEvent, UInt_t nRun, 
                         Int_t sourceId, UInt_t sequence, Double_t timestamp,
                         const TString& str) {
    m_eventId = nEvent;
    m_runId = nRun;
    m_sourceId = sourceId;
    m_sequence = sequence;
    m_timeStamp = timestamp;
    m_sourceName = str;
    return;
}


void McEvent::Clear(Option_t *option) {
    

    m_eventId = 0;
    m_runId = 0;
    m_sourceName = sourceUnknownStr;
    
    if (option) {
        TString optStr(option);
        if (optStr.CompareTo("ALL")==0) {
            if (m_particleCol)       m_particleCol->Clear();
            if (m_positionHitCol)    m_positionHitCol->Clear();
            if (m_integratingHitCol) m_integratingHitCol->Clear();
            if (m_trajectoryCol)     m_trajectoryCol->Clear();
            return;
        }
    }

    // Tell the MC object manager to reset its iterators
    McObjectManager::getPointer()->Delete();

    if (m_particleCol) {
        m_particleCol->Clear();
    }

    if (m_positionHitCol) {
        m_positionHitCol->Clear();
   }

    if (m_tkrStripCol)
    {
        m_tkrStripCol->SetOwner();
        m_tkrStripCol->Clear();
    }

    if (m_integratingHitCol) {
        m_integratingHitCol->Clear();
    }

    if (m_trajectoryCol) 
    {
        m_trajectoryCol->Clear();
    }
}

void McEvent::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    std::cout.precision(2);
    std::cout << "Run: " << m_runId << " Event: " << m_eventId << std::endl;
    std::cout << "SourceName: " << m_sourceName << std::endl;
    std::cout << m_particleCol->GetEntries() << " McParticles" << std::endl;
    std::cout << m_positionHitCol->GetEntries() << " McPositionHits" << std::endl;
    std::cout << m_integratingHitCol->GetEntries() 
        << " McIntegratingHits" << std::endl;
    std::cout << m_trajectoryCol->GetEntries() << " McTrajectories" << std::endl;
}


void McEvent::addMcParticle(McParticle *part) {   
    m_particleCol->Add(part);
    return;
}

McParticle* McEvent::getMcParticle(UInt_t index) const {
    return ((McParticle*)m_particleCol->At(index));
}

void McEvent::addMcPositionHit(McPositionHit *hit) {
    m_positionHitCol->Add(hit);
    return;
}

McPositionHit* McEvent::getMcPositionHit(UInt_t index) const {

    return ((McPositionHit*)m_positionHitCol->At(index));
}

void McEvent::addMcTkrStrip(McTkrStrip *strip) {
    m_tkrStripCol->Add(strip);
    return;
}

McTkrStrip* McEvent::getMcTkrStrip(UInt_t index) const {

    return ((McTkrStrip*)m_tkrStripCol->At(index));
}

void McEvent::addMcIntegratingHit(McIntegratingHit *hit) {
    m_integratingHitCol->Add(hit);
    return;
}

McIntegratingHit* McEvent::getMcIntegratingHit(UInt_t index) const {
    return ((McIntegratingHit*)m_integratingHitCol->At(index));
}

void McEvent::addMcTrajectory(McTrajectory *hit) {
    m_trajectoryCol->Add(hit);
    return;
}

McTrajectory* McEvent::getMcTrajectory(UInt_t index) const {

    return ((McTrajectory*)m_trajectoryCol->At(index));
}


//======================================================
// For Unit Tests
//======================================================

// DC: I CANNOT MAKE THE USUAL Fake() and CompareInRange(),
// because McEvent is designed as a singleton.

const UInt_t NUM_ENTRIES = 5 ;
const Int_t SOURCE_ID = 7 ;
const UInt_t SEQUENCE = 4 ;
    
void McEvent::Fake( Int_t ievent, Int_t irun, Float_t randNum ) {


    Clear() ;
    
    //DC : it is strange that in McEvent ievent and irun
    //are unsigned, and not in the usual Fake. I should
    //have a look at all kinds of events, and either
    //change McEvent or Fake
    initialize(ievent,irun,SOURCE_ID,SEQUENCE,ievent*1.0, sourceUnknownStr) ;
        
    UInt_t ientry ;
    for ( ientry=0 ; ientry<NUM_ENTRIES ; ++ientry ) {
            
        McParticle * mcPart = new McParticle();
        mcPart->Fake(ievent,ientry,randNum) ;
        addMcParticle(mcPart) ;
            
        McPositionHit * posHit = new McPositionHit() ;
        posHit->Fake(ievent,ientry,randNum) ;
        addMcPositionHit(posHit) ;
            
        McIntegratingHit *intHit = new McIntegratingHit();
        intHit->Fake(ievent,ientry,randNum) ;
        intHit->addEnergyItem(1.5,mcPart,mcPart->getFinalPosition()) ;
        addMcIntegratingHit(intHit) ;
        
        McTrajectory *trajectory = new McTrajectory() ;
        trajectory->Fake(ievent,ientry,randNum) ;
        addMcTrajectory(trajectory) ;
    }
}

Bool_t McEvent::CompareToFake( Int_t iEvent, Int_t iRun, Float_t randNum ) {

    bool result = true ;
        
    result = rootdatautil::CompareInRange(getEventId(),(UInt_t)iEvent,"EventId") && result ;
    result = rootdatautil::CompareInRange(getRunId(),(UInt_t)iRun,"RunId") && result ;
    result = rootdatautil::CompareInRange(getSourceId(),SOURCE_ID,"SourceId") && result ;
    result = rootdatautil::CompareInRange(getSequence(),SEQUENCE,"Sequence") && result ;
    result = rootdatautil::CompareInRange(getTimeStamp(),iEvent*1.0,"TimeStamp") && result ;
    result = rootdatautil::CompareInRange(getSourceName(),sourceUnknownStr.Data(),"SourceName") && result ;
        
    TObjArray refParticleCol ;
    TObjArray refPositionHitCol ;
    TObjArray refIntHitCol ;
    TObjArray refTrajectoryCol ;

    UInt_t iEntry ;
    for ( iEntry=0 ; iEntry<NUM_ENTRIES ; ++iEntry ) {
        
        McParticle * particle = getMcParticle(iEntry) ;
        
        McParticle * mcPart = new McParticle();
        mcPart->Fake(iEvent,iEntry,randNum) ;
        refParticleCol.Add(mcPart) ;
        
        McPositionHit * posHit = new McPositionHit() ;
        posHit->Fake(iEvent,iEntry,randNum) ;
        refPositionHitCol.Add(posHit) ;
        
        McIntegratingHit *intHit = new McIntegratingHit();
        intHit->Fake(iEvent,iEntry,randNum) ;
        intHit->addEnergyItem(1.5,particle,particle->getFinalPosition()) ;
        refIntHitCol.Add(intHit) ;
        
        McTrajectory *trajectory = new McTrajectory() ;
        trajectory->Fake(iEvent,iEntry,randNum) ;
        refTrajectoryCol.Add(trajectory) ;
    
    }
    
    result = rootdatautil::TObjArrayCompareInRange<McParticle>(getMcParticleCol(),&refParticleCol,"McParticle") && result ;
    result = rootdatautil::TObjArrayCompareInRange<McPositionHit>(getMcPositionHitCol(),&refPositionHitCol,"McPositionHit") && result ;
    result = rootdatautil::TObjArrayCompareInRange<McIntegratingHit>(getMcIntegratingHitCol(),&refIntHitCol,"McIntegratingHit") && result ;
    result = rootdatautil::TObjArrayCompareInRange<McTrajectory>(getMcTrajectoryCol(),&refTrajectoryCol,"McTrajectory") && result ;

    std::string name = "McEvent" ;
        
    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}

//#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)
//#define COMPARE_TOBJ_ARRAY_IN_RANGE(T) rootdatautil::TObjArrayCompareInRange<T>(get ## T ## Col(),ref.get ## T ## Col(),#T)
//
//Bool_t McEvent::CompareInRange( const McEvent & ref, const std::string & name ) {
//
//    bool result = true ;
//
//    result = COMPARE_IN_RANGE(EventId) && result ;
//    result = COMPARE_IN_RANGE(RunId) && result ;
//    result = COMPARE_IN_RANGE(SourceId) && result ;
//    result = COMPARE_IN_RANGE(Sequence) && result ;
//    result = COMPARE_IN_RANGE(TimeStamp) && result ;
//
//    result = COMPARE_TOBJ_ARRAY_IN_RANGE(McParticle) && result ;
//    result = COMPARE_TOBJ_ARRAY_IN_RANGE(McPositionHit) && result ;
//    result = COMPARE_TOBJ_ARRAY_IN_RANGE(McIntegratingHit) && result ;
//    result = COMPARE_TOBJ_ARRAY_IN_RANGE(McTrajectory) && result ;
//
//    if (!result) {
//        if ( name == "" ) {
//            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
//        }
//        else {
//            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
//        }
//    }
//    return result ;
//
//}

