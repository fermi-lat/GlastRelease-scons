#include "digiRootData/CalHeader.h"


ClassImp(CalHeader)


CalHeader::CalHeader() {
    m_eventId = 0;
    m_timerWord = 0;
    m_TREQ_VETO_status = 0;
    m_deadtime_and_cause = 0;
}

CalHeader::~CalHeader() {
}

void CalHeader::Clean(Option_t *option) {
    m_eventId = 0;
    m_timerWord = 0;
    m_TREQ_VETO_status = 0;
    m_deadtime_and_cause = 0;
}