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
	m_mcDaugterParticles = new TObjArray();
    m_mcDaugterParticles->SetOwner();
}
//________________________________________________________________________
RMcVertex::~RMcVertex() {
    if (m_pPosI)
		delete m_pPosI;
	
    if (m_pPosF)
		delete m_pPosF;
	
    if (m_pMomI)
		delete m_pMomI;
	
    if (m_pMomF)
		delete m_pMomF;

	if (m_mcDaugterParticles)
		delete m_mcDaugterParticles;
}
//________________________________________________________________________

void RMcVertex::addDaughterMcParticle( RMcParticle* particle )
{
	if (particle) {
        m_mcDaugterParticles->Add(particle);
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
    return m_mcDaugterParticles;
}