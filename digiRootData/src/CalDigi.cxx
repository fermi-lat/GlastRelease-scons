#include "digiRootData/CalDigi.h"
#include "digiRootData/CalHit.h"
#include "digiRootData/LogID.h"


ClassImp(CalDigi)

TClonesArray *CalDigi::m_staticCalArray = 0;


CalDigi::CalDigi() {
    if (!m_staticCalArray) m_staticCalArray = new TClonesArray("CalHit", 80);
    m_logs = m_staticCalArray;
    m_numLogs = -1;
    m_eventId = 0;
    m_timerWord = 0;
    m_TREQ_VETO_status = 0;
    m_deadTime = 0;
    m_deadTimeCause = 0;
}

CalDigi::~CalDigi() {
    // Destructor
    m_logs->Delete();
    if (m_logs == m_staticCalArray) m_staticCalArray = 0;
    delete m_logs;
}

CalHit* CalDigi::Add(LogID *id) {
    ++m_numLogs;
    TClonesArray &logs = *m_logs;
    new(logs[m_numLogs]) CalHit(id);
    return ((CalHit*)(logs[m_numLogs]));
}

void CalDigi::Clean(Option_t *option) {
    // Call Delete rather than Clear
    // since the elements in this array
    // have allocated memory
    m_logs->Delete(option);
}