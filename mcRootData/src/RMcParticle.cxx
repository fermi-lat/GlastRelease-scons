///////////////////////////////////////////////////////////////////////////
//                                                                       
// The RMcParticle class provides an interface to a particle for use in
// GLAST Monte Carlo applications
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "mcRootData/RMcParticle.h"
#include "mcRootData/RMcVertex.h"

ClassImp(RMcParticle)
//________________________________________________________________________
/*! Default Constructor
 *  Reserved for use by ROOT I/O.  Do not use this constructor, as no
 *  initialization is performed.
 */
RMcParticle::RMcParticle() :
    m_particleId(0),  m_mcVertex(0) 
{
}
//________________________________________________________________________
RMcParticle::RMcParticle(Int_t id, RMcVertex *pVertex) :
  m_particleId(id), m_mcVertex(pVertex)    
{
      m_primary =0;
}
//________________________________________________________________________
RMcParticle::~RMcParticle() {
  if (m_mcVertex)
    delete m_mcVertex;

}

RMcVertex* RMcParticle::mcVertex() {
	return m_mcVertex;
}


void RMcParticle::setMcVertex(RMcVertex* value){
	m_mcVertex = value;
}

void RMcParticle::setParticleId(Int_t id)
{
    m_particleId = id;
}

Int_t RMcParticle::particleId()
{
    return m_particleId;
}

void RMcParticle::setPrimaryParticleFlag(Bool_t isPrimary)
{
    m_primary = isPrimary;
}

Bool_t RMcParticle::primaryParticle()
{
    return m_primary;
}