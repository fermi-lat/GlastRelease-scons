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
    m_nEvent(0), m_nRun(0), m_nVertCount(0), m_pVertList(0)
{
}
//________________________________________________________________________
RMcEvent::~RMcEvent() {
    Clean();
}
//________________________________________________________________________
void RMcEvent::Clean() {
    if (m_pVertList)
        delete m_pVertList;

    m_pVertList = 0;

    m_nEvent = 0;
    m_nRun = 0;
    m_nVertCount = 0;
}
//________________________________________________________________________
void RMcEvent::Create() {
    m_pVertList = new TObjArray();
    m_pVertList->SetOwner();

    m_nEvent = 0;
    m_nRun = 0;
    m_nVertCount = 0;
}
//________________________________________________________________________
/*! Returns requested particle if it exists.  Else returns NULL.
 */
RMcVertex *RMcEvent::getVert(Int_t index) const {
    if ((index >= 0) && (index < getVertCount()))
        return (RMcVertex*)m_pVertList->At(index);
    else
        return 0;
}
//________________________________________________________________________
/*! Adds particle (if non-NULL) to list of particles
 *  This function, and ALL functions which modify the particle list
 *  MUST also update the m_nVertCount member.
 */
Int_t RMcEvent::addVert(RMcVertex *pVert) {
    if (pVert) {
        m_pVertList->Add(pVert);
        return (m_nVertCount++);
    } else
        return -1;
}