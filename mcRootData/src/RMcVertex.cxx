///////////////////////////////////////////////////////////////////////////
//                                                                       
// The RMcVertex class provides a position and momentum information for
// inclusion in the McParticle class.
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "mcRootData/RMcVertex.h"
#include "mcRootData/RMcParticle.h"

ClassImp(RMcVertex)
//________________________________________________________________________
/*! Default Constructor
*  Reserved for use by ROOT I/O.  Do not use this constructor, as no
*  initialization is performed.
*/
RMcVertex::RMcVertex() :
m_pPosI(0), m_pPosF(0), m_pMomI(0), m_pMomF(0)
{
}
//________________________________________________________________________
RMcVertex::RMcVertex(TLorentzVector *pPosI, TLorentzVector *pPosF, 
					 TLorentzVector *pMomI, TLorentzVector *pMomF) :
m_pPosI(pPosI), m_pPosF(pPosF), m_pMomI(pMomI), m_pMomF(pMomF)
{
	m_mcDaughterParticles = new TObjArray();
}
//________________________________________________________________________
RMcVertex::~RMcVertex() {
  
	if (m_pPosI){
		delete m_pPosI;
		m_pPosI = 0;
	}
	
    if (m_pPosF){
		delete m_pPosF;
		m_pPosF =0;
	}
	
    if (m_pMomI){
		delete m_pMomI;
		m_pMomI = 0;
	}
	
    if (m_pMomF){
		delete m_pMomF;
		m_pMomF = 0;
	}

	if (m_mcDaughterParticles){
		delete m_mcDaughterParticles;
		m_mcDaughterParticles = 0;
	}
	
}
//________________________________________________________________________

void RMcVertex::addDaughterMcParticle( RMcParticle* particle )
{
	if (particle) {
        m_mcDaughterParticles->Add(particle);
    }
}


void  RMcVertex::setMotherMcParticle(RMcParticle* mother)
{
	m_motherParticle = mother;
}

RMcParticle* RMcVertex::mcParticle()
{
    return m_mainParticle;
}

void RMcVertex::setMcParticle(RMcParticle* value)
{
    m_mainParticle = value;
}

RMcParticle* RMcVertex::motherMcParticle()
{
    return m_motherParticle;
}

const TObjArray* RMcVertex::daughterMcParticles() const
{
    return m_mcDaughterParticles;
}