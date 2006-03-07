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
#include <iostream>

// Allocate the TObjArrays just once
TObjArray *McEvent::m_staticParticleCol = 0;
TObjArray *McEvent::m_staticPositionHitCol = 0;
TObjArray *McEvent::m_staticIntegratingHitCol = 0;
TObjArray *McEvent::m_staticTrajectoryCol = 0;

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
    
    if (!m_staticTrajectoryCol) m_staticTrajectoryCol = new TObjArray();
    m_trajectoryCol = m_staticTrajectoryCol;
}

McEvent::~McEvent() {
    
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

    if (m_trajectoryCol == m_staticTrajectoryCol) m_staticTrajectoryCol = 0;
    m_trajectoryCol->Delete();
    delete m_trajectoryCol;
    m_trajectoryCol = 0;

    Clear("ALL");
}

void McEvent::initialize(UInt_t nEvent, UInt_t nRun, 
                         Int_t sourceId, UInt_t sequence, Double_t timestamp) {
    m_eventId = nEvent;
    m_runId = nRun;
    m_sourceId = sourceId;
    m_sequence = sequence;
    m_timeStamp = timestamp;
    return;
}


void McEvent::Clear(Option_t *option) {
    

    m_eventId = 0;
    m_runId = 0;
    const Int_t ndpos = 150000;
    const Int_t nd = 150000;
    const Int_t ndtraj = 150000;
    static Int_t limit=100;
    static Int_t indpart=0;
    static Int_t indpos=0;
    static Int_t indint=0;
    static Int_t indtraj=0;
    static McParticle* keeppart[nd];
    static McPositionHit* keeppos[ndpos];
    static McIntegratingHit* keepint[nd];
    static McTrajectory* keeptraj[nd];
    

    if ( (option) && (option[0] == 'A') ) {
        Int_t j;
        for (j=0;j<indpart;j++) delete keeppart[j];
        indpart = 0;
        for (j=0;j<indpos;j++) delete keeppos[j];
        indpos = 0;
	for (j=0;j<indint;j++) delete keepint[j];
	indint = 0;
        for (j=0;j<indtraj;j++) delete keeptraj[j];
        indtraj = 0;
        if (m_particleCol) m_particleCol->Delete();
        if (m_positionHitCol) m_positionHitCol->Delete();
        if (m_integratingHitCol) m_integratingHitCol->Delete();
        if (m_trajectoryCol) m_trajectoryCol->Delete();
        return;
    }


    if (m_particleCol) {
      //      m_particleCol->Delete();
      Int_t n = m_particleCol->GetEntries();
      if (n>limit) {
        //std::cout <<"!!!Warning: particle nr entries more than limit!!!Limit was increased "<<n<<std::endl;
	limit=n+10;
        if (limit > nd)
          std::cout << "!!!Warning:  limit for mcParticles is greater than " << nd << std::endl;
	for (Int_t j=0;j<indpart;j++) delete keeppart[j];
	indpart = 0;
      }
      for (Int_t i=0;i<n;i++) keeppart[indpart+i] = (McParticle*)m_particleCol->At(i);
      indpart += n;
      if (indpart > nd-limit) {
	for (Int_t j=0;j<indpart;j++) delete keeppart[j];
	indpart = 0;
      }
      m_particleCol->Clear();
    }

    if (m_positionHitCol) {
      //      m_positionHitCol->Delete();
      Int_t n = m_positionHitCol->GetEntries();
      if (n>limit) {
        //std::cout <<"!!!Warning: positionhit nr entries more than limit!!!Limit was increased "<<n<<std::endl;
	limit=n+10;
        if (limit > ndpos)
          std::cout << "!!!Warning:  limit for posHits is greater than " << ndpos << std::endl;
        for (Int_t j=0;j<indpos;j++) delete keeppos[j];
	indpos = 0;
      }

      for (Int_t i=0;i<n;i++) keeppos[indpos+i] = (McPositionHit*)m_positionHitCol->At(i);
      indpos += n;
      if (indpos > ndpos-limit) {
        for (Int_t j=0;j<indpos;j++) delete keeppos[j];
	indpos = 0;
      }
      m_positionHitCol->Clear();
      m_positionHitCol->Expand(10000);  //temporary
   }

    if (m_integratingHitCol) {
      //      m_integratingHitCol->Delete();
      Int_t n = m_integratingHitCol->GetEntries();
      if (n>limit) {
        //std::cout <<"!!!Warning: integratinghit nr entries more than limit!!!Limit was increased "<<n<<std::endl;
	limit=n+10;
        if (limit > nd)
          std::cout << "!!!Warning:  limit for intHits is greater than " << nd << std::endl;
	for (Int_t j=0;j<indint;j++) delete keepint[j];
	indint = 0;
      }
      for (Int_t i=0;i<n;i++) keepint[indint+i] = (McIntegratingHit*)m_integratingHitCol->At(i);
      indint += n;
      if (indint > nd-limit) {
	for (Int_t j=0;j<indint;j++) delete keepint[j];
	indint = 0;
      }
      m_integratingHitCol->Clear();
    }

    if (m_trajectoryCol) 
    {
        Int_t n = m_trajectoryCol->GetEntries();
        if (n>limit) 
        {
            //std::cout <<"!!!Warning: positionhit nr entries more than limit!!!Limit was increased "<<n<<std::endl;
	        limit=n+10;
            if (limit > ndtraj)
                std::cout << "!!!Warning:  limit for McTrajectory is greater than " << ndtraj << std::endl;
            for (Int_t j=0;j<indtraj;j++) delete keeptraj[j];
	        indtraj = 0;
        }

        for (Int_t i=0;i<n;i++) keeptraj[indtraj+i] = (McTrajectory*)m_trajectoryCol->At(i);
        indtraj += n;
        if (indtraj > ndtraj-limit) 
        {
            for (Int_t j=0;j<indtraj;j++) delete keeptraj[j];
	        indtraj = 0;
        }
        m_trajectoryCol->Clear();
        m_trajectoryCol->Expand(10000);  //temporary
    }
}

void McEvent::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    std::cout.precision(2);
    std::cout << "Run: " << m_runId << " Event: " << m_eventId << std::endl;
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
