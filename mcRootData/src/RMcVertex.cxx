///////////////////////////////////////////////////////////////////////////
//                                                                       
// The RMcVertex class provides a position and momentum information for
// inclusion in the McParticle class.
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "mcRootData/RMcVertex.h"

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
}
//________________________________________________________________________
