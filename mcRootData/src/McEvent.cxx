///////////////////////////////////////////////////////////////////////////
//                                                                       
// The McEvent class contains general information about a Monte Carlo
// Event, (Run/Event numbers, etc) as well as a list of McParticle
// instances generated during the simulated event.
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "mcRootData/McEvent.h"

ClassImp(McEvent)
//________________________________________________________________________
/*! Default Constructor
 *  NORMALLY Reserved for use by ROOT I/O.  For this class, all initialization
 *  is performed by calling the Clean / Create methods and filling the data
 *  members using the setXXXX functions provided.  Default constructor is,
 *  therefore, used for instantiating new instances.
 */
McEvent::McEvent() : 
    m_nEvent(0), m_nRun(0), m_nPartCount(0), m_pPartList(0)
{
}
//________________________________________________________________________
McEvent::~McEvent() {
    Clean();
}
//________________________________________________________________________
void McEvent::Clean() {
    if (m_pPartList)
        delete m_pPartList;

    m_pPartList = 0;

    m_nEvent = 0;
    m_nRun = 0;
    m_nPartCount = 0;
}
//________________________________________________________________________
void McEvent::Create() {
    m_pPartList = new TObjArray();
    m_pPartList->SetOwner();

    m_nEvent = 0;
    m_nRun = 0;
    m_nPartCount = 0;
}
//________________________________________________________________________
/*! Returns requested particle if it exists.  Else returns NULL.
 */
McParticle *McEvent::getPart(Int_t index) const {
    if ((index >= 0) && (index < getPartCount()))
        return (McParticle*)m_pPartList->At(index);
    else
        return 0;
}
//________________________________________________________________________
/*! Adds particle (if non-NULL) to list of particles
 *  This function, and ALL functions which modify the particle list
 *  MUST also update the m_nPartCount member.
 */
Int_t McEvent::addPart(McParticle *pPart) {
    if (pPart) {
        m_pPartList->AddLast(pPart);
        return (m_nPartCount++);
    } else
        return -1;
}