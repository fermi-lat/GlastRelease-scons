///////////////////////////////////////////////////////////////////////////
//                                                                       
// The McParticle class provides an interface to a particle for use in
// GLAST Monte Carlo applications
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "mcRootData/McParticle.h"

ClassImp(McParticle)
//________________________________________________________________________
/*! Default Constructor
 *  Reserved for use by ROOT I/O.  Do not use this constructor, as no
 *  initialization is performed.
 */
McParticle::McParticle() :
    m_nPDGId(0), m_pParent(0), m_pVert(0) 
{
}
//________________________________________________________________________
McParticle::McParticle(Int_t id, McParticle *pParent, McVertex *pVertex) :
  m_nPDGId(id), m_pParent(pParent), m_pVert(pVertex)    
{
    m_pChildList = new TObjArray();
    m_pChildList->SetOwner();
}
//________________________________________________________________________
McParticle::~McParticle() {
  if (m_pVert)
    delete m_pVert;

  if (m_pChildList)
      delete m_pChildList;
}
//________________________________________________________________________
/*! Returns requested child if it exists.  Else returns NULL.
 */
McParticle *McParticle::getChild(Int_t nIndex) const {
    McParticle *pChild = 0;
    if ((nIndex >= 0) && (nIndex < m_nChildCount))
        pChild = (McParticle*)m_pChildList->At(nIndex);
    return pChild;
}
//________________________________________________________________________
/*! Adds child (if non-NULL) to list of children
 *  This function, and ALL functions which modify the child list
 *  MUST also update the m_nChildCount member.
 */
void McParticle::addChild(McParticle *pChild) {
    if (pChild) {
        m_pChildList->Add(pChild);
        m_nChildCount++;
    }
}
//________________________________________________________________________