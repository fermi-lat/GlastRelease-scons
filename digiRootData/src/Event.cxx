///////////////////////////////////////////////////////////////////////////
//
// The Event class is the top level class for the ASCII to ROOT
// conversion. It contains all the information about a single event.
//
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/Event.h"

ClassImp(Event)

CalDigi *Event::m_staticCal=0;
AcdDigi *Event::m_staticAcd=0;
TkrDigi *Event::m_staticTkr=0;
L1T *Event::m_staticL1T=0;

//__________________________________________________________________________
Event::Event() {
    // Default constructor.
    // Assign default values to members
    
   // if (!m_staticAcd) m_staticAcd = new AcdDigi();
   // m_ACD = m_staticAcd;
    
    //if (!m_staticCal) m_staticCal = new CalDigi();
   // m_CAL = m_staticCal;
    
    //if (!m_staticTkr) m_staticTkr = new TkrDigi();
    //m_TKR = m_staticTkr;
    
    //if (!m_staticL1T) m_staticL1T = new L1T();
    //m_L1Trigger = m_staticL1T;
    
    m_ACD = 0;
    m_CAL = 0;
    m_TKR = 0;
    m_L1Trigger = 0;

  //  m_Tagger = 0;  
    m_eventId = 0;
    m_run = 0;
}
//_________________________________________________________________________
Event::~Event() {
  // Destructor
    Clean();
    
    //if (m_ACD == m_staticAcd) m_staticAcd = 0;
    if (m_ACD)
    delete m_ACD;
    //if (m_CAL == m_staticCal) m_staticCal = 0;
    if (m_CAL)
    delete m_CAL;
    //if (m_TKR == m_staticTkr) m_staticTkr = 0;
    if (m_TKR)
    delete m_TKR;
    //if (m_L1Trigger == m_staticL1T) m_staticL1T = 0;
    if (m_L1Trigger)
    delete m_L1Trigger;

}
//_________________________________________________________________________
void Event::Clean(Option_t *option) {
    if (m_ACD) {
        m_ACD->Clean();
        delete m_ACD;
        m_ACD = 0;
    }
    if (m_CAL) {
        m_CAL->Clean();
        delete m_CAL;
        m_CAL = 0;
    }
    if (m_TKR) {
        m_TKR->Clean();
        delete m_TKR;
        m_TKR = 0;
    }
    if (m_L1Trigger) {
        delete m_L1Trigger;
        m_L1Trigger = 0;
    }
}


