#include "reconRootData/ReconEvent.h"
#include <iostream>

ClassImp(ReconEvent)

ReconEvent::ReconEvent() {
    Clear();
}

ReconEvent::~ReconEvent() {
}

void ReconEvent::initialize(UInt_t eventId, UInt_t runId) {
    m_eventId = eventId;
    m_runId = runId;
}

void ReconEvent::Clear(Option_t *option) {
    m_eventId = 0;
    m_runId = 0;
}

void ReconEvent::Print(Option_t *option) {
    using namespace std;
    cout << "(Run, Event): (" << m_runId << ", " << m_eventId << ")" << endl;
}
