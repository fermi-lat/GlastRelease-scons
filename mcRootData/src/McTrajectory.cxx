// File and Version Information:
//      $Header$
//
// Description:
//  The McTrajectory class contains the step by step position information for
//  particles tracked through GLAST during the simulation stage. 
//                                                
// Author:  Tracy Usher                    

#include "mcRootData/McTrajectory.h"
#include "Riostream.h"

ClassImp(McTrajectory)

McTrajectory::McTrajectory() : m_mcParticle(0), m_charge(0)
{
    m_pointCol = new TObjArray();

    m_pointCol->Clear();
}

McTrajectory::~McTrajectory() 
{    
    m_pointCol->Delete();
    delete m_pointCol;
    m_pointCol = 0;

	Clear();
}

void McTrajectory::initialize(McParticle* mcPart, UInt_t charge) 
{
    m_mcParticle = mcPart;
    m_charge     = charge;

    return;
}


void McTrajectory::Clear(Option_t *option) 
{
    m_mcParticle = 0;
    m_charge     = 0;

    if (m_pointCol) m_pointCol->Clear();

    return;
}

void McTrajectory::Print(Option_t *option) const 
{
/*
    using namespace std;
    TObject::Print(option);
    std::cout.precision(2);
    std::cout << "Run: " << m_runId << " Event: " << m_eventId << std::endl;
    std::cout << m_particleCol->GetEntries() << " McParticles" << std::endl;
    std::cout << m_positionHitCol->GetEntries() << " McPositionHits" << std::endl;
    std::cout << m_integratingHitCol->GetEntries() 
        << " McIntegratingHits" << std::endl;
*/
}

