#include "digiRootData/CalHeader.h"


ClassImp(CalHeader)


CalHeader::CalHeader() {
    m_eventId = 0;
    m_timerWord = 0;
    m_TREQ_VETO_status = 0;
    m_deadTime = 0;
    m_deadTimeCause = 0;
}

CalHeader::~CalHeader() {
}

void CalHeader::setHeader(UInt_t id, UInt_t time, UInt_t stat, 
                          UInt_t deadTime, UInt_t cause) {
    m_eventId = id;
    m_timerWord = time;
    m_TREQ_VETO_status = stat;
    m_deadTime = deadTime;
    m_deadTimeCause = cause;
}

void CalHeader::Clean(Option_t *option) {
    m_eventId = 0;
    m_timerWord = 0;
    m_TREQ_VETO_status = 0;
    m_deadTime = 0;
    m_deadTimeCause = 0;
}