// File and Version Information:
//      $Header$
//
// Description:
// The McEvent class contains general information about a Monte Carlo
// Event, (Run/Event numbers, etc) as well as a list of McParticle
// instances generated during the simulated event.
//                                                
// Authors:  Dan Flath and Heather Kelly                       

#include "mcRootData/McEvent.h"

// Allocate the TObjArrays just once
TObjArray *McEvent::m_staticParticleCol = 0;
TObjArray *McEvent::m_staticPositionHitCol = 0;
TObjArray *McEvent::m_staticIntegratingHitCol = 0;

ClassImp(McEvent)

McEvent::McEvent() : 
m_eventId(0), m_runId(0)
{
    if (!m_staticParticleCol) m_staticParticleCol = new TObjArray();
    m_particleCol = m_staticParticleCol;
    
    if (!m_staticPositionHitCol) m_staticPositionHitCol = new TObjArray();
    m_positionHitCol = m_staticPositionHitCol;
    
    if (!m_staticIntegratingHitCol) m_staticIntegratingHitCol = new TObjArray();
    m_integratingHitCol = m_staticIntegratingHitCol;
}

McEvent::~McEvent() {
    Clear();
    
    if (m_particleCol == m_staticParticleCol) m_staticParticleCol = 0;
    m_particleCol->Delete();
    delete m_particleCol;
    m_particleCol = 0;

    if (m_positionHitCol == m_staticPositionHitCol) m_staticPositionHitCol = 0;
    m_positionHitCol->Delete();
    delete m_positionHitCol;
    m_positionHitCol = 0;

    if(m_integratingHitCol == m_staticIntegratingHitCol) m_staticIntegratingHitCol = 0;
    m_integratingHitCol->Delete();
    delete m_integratingHitCol;
    m_integratingHitCol = 0;
    
}

void McEvent::initialize(Int_t nEvent, Int_t nRun) {
    m_eventId = nEvent;
    m_runId = nRun;
    return;
}


void McEvent::Clear(Option_t *option) {
    
    m_eventId = 0;
    m_runId = 0;
    
    m_particleCol->Clear();
    m_positionHitCol->Clear();
    m_integratingHitCol->Clear();
}


void McEvent::addMcParticle(McParticle *part) {   
    m_particleCol->Add(part);
    return;
}

McParticle* McEvent::getMcParticle(Int_t index) const {
    return ((McParticle*)m_particleCol->At(index));
}

void McEvent::addMcPositionHit(McPositionHit *hit) {
    m_positionHitCol->Add(hit);
    return;
}

McPositionHit* McEvent::getMcPositionHit(Int_t index) const {

    return ((McPositionHit*)m_positionHitCol->At(index));
}

void McEvent::addMcIntegratingHit(McIntegratingHit *hit) {
    m_integratingHitCol->Add(hit);
    return;
}

McIntegratingHit* McEvent::getMcIntegratingHit(Int_t index) const {
    return ((McIntegratingHit*)m_integratingHitCol->At(index));
}