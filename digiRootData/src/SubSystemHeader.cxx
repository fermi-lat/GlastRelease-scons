#include "digiRootData/SubSystemHeader.h"


ClassImp(SubSystemHeader)

SubSystemHeader::SubSystemHeader() {
    m_eventId = 0;
    m_timerWord = 0;
    m_TREQ_VETO_status = 0;
    m_deadtime_and_cause = 0;
}

SubSystemHeader::~SubSystemHeader() {
}

void SubSystemHeader::Clean(Option_t *option) {
    m_eventId = 0;
    m_timerWord = 0;
    m_TREQ_VETO_status = 0;
    m_deadtime_and_cause = 0;
}