#include "reconRootData/TkrGamma.h"

ClassImp(TkrGamma)

TkrGamma::TkrGamma() {
    // default constructor
    m_tracks = 0;
    m_locator = 0;
}

TkrGamma::TkrGamma(Int_t dummy) {
    m_tracks = 0;
    m_locator = 0;
    Create();
}


TkrGamma::~TkrGamma() {
    Clean();
    if (m_locator) {
        delete m_locator;
        m_locator = 0;
    }
    if (m_tracks) {
        delete m_tracks;
        m_tracks = 0;
    }
}

void TkrGamma::Clean() {
    // Does not de-allocate the tracks, since
    // probably pointers into TkrRecon's list of
    // tracks
    // But do free the TObjArray
    if (m_tracks) m_tracks->Clear();

    // but we do free memory associated with 
    // TkrLocators for this gamma
    if (m_locator) {
        int nEntries = m_locator->GetEntries();
        for (int i=0; i<nEntries; i++)
            delete m_locator->At(i);
        m_locator->Clear();
        //delete m_locator;
        //m_locator = 0;
    }
}

void TkrGamma::Create() {
    if (!m_locator) m_locator = new TObjArray();
    if (!m_tracks) m_tracks = new TObjArray();
}
