///////////////////////////////////////////////////////////////////////////
//                                                                       
// The RMcParticle class provides an interface to a particle for use in
// GLAST Monte Carlo applications
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "mcRootData/RMcParticle.h"

ClassImp(RMcParticle)
//________________________________________________________________________
/*! Default Constructor
 *  Reserved for use by ROOT I/O.  Do not use this constructor, as no
 *  initialization is performed.
 */
RMcParticle::RMcParticle() :
    m_nPDGId(0), m_pParent(0), m_pVert(0) 
{
}
//________________________________________________________________________
RMcParticle::RMcParticle(Int_t id, RMcParticle *pParent, RMcVertex *pVertex) :
  m_nPDGId(id), m_pParent(pParent), m_pVert(pVertex)    
{
    m_pChildList = new TObjArray();
    m_pChildList->SetOwner();
}
//________________________________________________________________________
RMcParticle::~RMcParticle() {
  if (m_pVert)
    delete m_pVert;

  if (m_pChildList)
      delete m_pChildList;
}
//________________________________________________________________________
/*! Returns requested child if it exists.  Else returns NULL.
 */
RMcParticle *RMcParticle::getChild(Int_t nIndex) const {
    RMcParticle *pChild = 0;
    if ((nIndex >= 0) && (nIndex < m_nChildCount))
        pChild = (RMcParticle*)m_pChildList->At(nIndex);
    return pChild;
}
//________________________________________________________________________
/*! Adds child (if non-NULL) to list of children
 *  This function, and ALL functions which modify the child list
 *  MUST also update the m_nChildCount member.
 */
void RMcParticle::addChild(RMcParticle *pChild) {
    if (pChild) {
        m_pChildList->Add(pChild);
        m_nChildCount++;
    }
}
//________________________________________________________________________