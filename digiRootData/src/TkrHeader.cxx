#include "digiRootData/TkrHeader.h"


ClassImp(TkrHeader)

TkrHeader::TkrHeader() {
    m_eventId = 0;
    m_timerWord = 0;
    m_TREQ_VETO_status = 0;
    m_deadTime = 0;
    m_deadTimeCause = 0;
}

TkrHeader::~TkrHeader() {
}


void TkrHeader::Clean(Option_t *option) {
    m_eventId = 0;
    m_timerWord = 0;
    m_TREQ_VETO_status = 0;
    m_deadTime = 0;
    m_deadTimeCause = 0;
}