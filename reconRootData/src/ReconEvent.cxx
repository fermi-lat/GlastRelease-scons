#include "reconRootData/Recon.h"
#include <iostream>

ClassImp(Recon)

Recon::Recon() {
    Clear();
}

Recon::~Recon() {
}

void Recon::initialize(UInt_t eventId, UInt_t runId) {
    m_eventId = eventId;
    m_runId = runId;
}

void Recon::Clear(Option_t *option) {
    m_eventId = 0;
    m_runId = 0;
}

void Recon::Print(Option_t *option) {
    using namespace std;
    cout << "(Run, Event): (" << m_runId << ", " << m_eventId << ")" << endl;
}
