#include "reconRootData/CalRecon.h"

ClassImp(CalRecon)

CalRecon::CalRecon() {
    m_logs = 0;
    m_clusters = 0;
}


CalRecon::~CalRecon() {
    Clean();
    if(m_clusters) {
        delete m_clusters;
        m_clusters = 0;
    }
    if (m_logs) {
        delete m_logs;
        m_logs = 0;
    }
}

void CalRecon::Clean() {

    if (m_clusters) {
        
        int nEntries = m_clusters->GetEntries();
        for (int i=0; i<nEntries; i++)
            delete m_clusters->At(i);
            
        m_clusters->Clear();
        //delete m_clusters;
        //m_clusters = 0;
    }

    if (m_logs) {
        int nEntries = m_logs->GetEntries();
        for (int i=0; i<nEntries; i++)
            delete m_logs->At(i);
        m_logs->Clear();
       // delete m_logs;
       // m_logs = 0;
    }
}

void CalRecon::Create() {
    if (!m_logs) m_logs = new TObjArray();
    if (!m_clusters) m_clusters = new TObjArray();
}
