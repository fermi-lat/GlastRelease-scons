#include "reconRootData/ReconEvent.h"
#include <iostream>

ClassImp(ReconEvent)

ReconEvent::ReconEvent() {
    m_tkr = 0;
    m_cal = 0;
    m_acd = 0;
    Clear();
}

ReconEvent::~ReconEvent() {
    Clear();

    if (m_cal) {
        delete m_cal;
        m_cal = 0;
    }

    if (m_tkr) {
        delete m_tkr;
        m_tkr = 0;
    }

    if (m_acd) {
        delete m_acd;
        m_acd = 0;
    }

}

void ReconEvent::initialize(UInt_t eventId, UInt_t runId, TkrRecon *tkr, CalRecon *cal, AcdRecon *acd) {
    m_eventId = eventId;
    m_runId = runId;
    m_tkr = tkr;
    m_cal = cal;
    m_acd = acd;
}

void ReconEvent::Clear(Option_t *option) {
    m_eventId = 0;
    m_runId = 0;
    if (m_cal) {
        delete m_cal;
        m_cal = 0;
    }
    if (m_tkr) {
        delete m_tkr;
        m_tkr = 0;
    }
    if (m_acd) {
        delete m_acd;
        m_acd = 0;
    }
}

void ReconEvent::Print(Option_t *option) const {
    using namespace std;
    cout << "(Run, Event): (" << m_runId << ", " << m_eventId << ")" << endl;
}
