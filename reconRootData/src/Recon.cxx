#include "reconRootData/Recon.h"

ClassImp(Recon)

// default constructor does not create
// any objects
Recon::Recon() {
   // m_tkr = 0;
   // m_cal = 0;
   // m_recFlags = 0;
//	m_name = "Recon";
}

Recon::~Recon() {
    Clean();
}

void Recon::Clean() {
    
    /*
    if (m_tkr) {
        delete m_tkr;
        m_tkr = 0;
    }
    
    if (m_cal) {
        delete m_cal;
        m_cal = 0;
    }

    if (m_recFlags) {
        delete m_recFlags;
        m_recFlags = 0;
    }
    */
}

void Recon::Create() {

    m_tkr.Create();
    m_cal.Create();
    /*
    if (!m_tkr) {
        m_tkr = new TkrRecon();
        m_tkr->Create();
    }

    if (!m_cal) {
        m_cal = new CalRecon();
        m_cal->Create();
    }

    if (!m_recFlags) {
        m_recFlags = new ReconHeader();
    }
    */
}
