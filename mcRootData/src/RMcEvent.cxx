///////////////////////////////////////////////////////////////////////////
//                                                                       
// The RMcEvent class contains general information about a Monte Carlo
// Event, (Run/Event numbers, etc) as well as a list of McParticle
// instances generated during the simulated event.
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "mcRootData/RMcEvent.h"

ClassImp(RMcEvent)
//________________________________________________________________________
/*! Default Constructor
 *  NORMALLY Reserved for use by ROOT I/O.  For this class, all initialization
 *  is performed by calling the Clean / Create methods and filling the data
 *  members using the setXXXX functions provided.  Default constructor is,
 *  therefore, used for instantiating new instances.
 */
RMcEvent::RMcEvent() : 
    m_nEvent(0), m_nRun(0), m_nPartCount(0), m_pPartList(0)
{
}
//________________________________________________________________________
RMcEvent::~RMcEvent() {
    Clean();
}
//________________________________________________________________________
void RMcEvent::Clean() {
    if (m_pPartList)
        delete m_pPartList;

    m_pPartList = 0;

    m_nEvent = 0;
    m_nRun = 0;
    m_nPartCount = 0;
}
//________________________________________________________________________
void RMcEvent::Create() {
    m_pPartList = new TObjArray();
    m_pPartList->SetOwner();

    m_nEvent = 0;
    m_nRun = 0;
    m_nPartCount = 0;
}
//________________________________________________________________________
/*! Returns requested particle if it exists.  Else returns NULL.
 */
RMcParticle *RMcEvent::getPart(Int_t index) const {
    if ((index >= 0) && (index < getPartCount()))
        return (RMcParticle*)m_pPartList->At(index);
    else
        return 0;
}
//________________________________________________________________________
/*! Adds particle (if non-NULL) to list of particles
 *  This function, and ALL functions which modify the particle list
 *  MUST also update the m_nPartCount member.
 */
Int_t RMcEvent::addPart(RMcParticle *pPart) {
    if (pPart) {
        m_pPartList->AddLast(pPart);
        return (m_nPartCount++);
    } else
        return -1;
}