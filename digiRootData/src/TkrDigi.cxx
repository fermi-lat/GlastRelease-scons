#include "digiRootData/TkrDigi.h"


ClassImp(TkrDigi)

TObjArray *TkrDigi::m_staticArray = 0;

TkrDigi::TkrDigi() {
    if (!m_staticArray) m_staticArray = new TObjArray();
    m_layers = m_staticArray;
    m_eventId = 0;
    m_timerWord = 0;
    m_TREQ_VETO_status = 0;
    m_deadTime = 0;
    m_deadTimeCause = 0;
}

TkrDigi::~TkrDigi() {
    // Destructor
    if (m_layers) m_layers->Delete();
 //   if (m_layers == m_staticArray) m_staticArray = 0;
 //   delete m_layers;
}

void TkrDigi::Add(TkrLayer *layer) {
    m_layers->Add(layer);
}

void TkrDigi::Clean(Option_t *option) {
    if (m_layers) m_layers->Delete(option);
}