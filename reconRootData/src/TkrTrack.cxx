#include "reconRootData/TkrTrack.h"

ClassImp(TkrTrack)

TkrTrack::TkrTrack() {
    // default constructor
    m_hits = 0;
    m_locator = 0;              // list of hits
    
    m_id = 0;              // track id
    m_type = NONE;            // type of the track (gamma, e+, e-, hadron ...)
    m_motherId = NONE;        // id of the mother particle (-1 if none)

    m_residual = 0.0f;         // residual  
    m_chi2 = 0.0f;             // chi**2 of the fit
    m_quality = 0.0f;          // quality measure of the track

    m_firstLayer = 0;         // layer where the particle starts

    m_energyInput = 0.0f;         // energy input from the cal
    m_energyDetermined = 0.0f;    // energy determined using the tracker + cal
}

TkrTrack::TkrTrack(UShort_t id) {    
    m_id = id;              // track id
    m_type = NONE;            // type of the track (gamma, e+, e-, hadron ...)
    m_motherId = NONE;        // id of the mother particle (-1 if none)

    m_residual = 0.0f;         // residual  
    m_chi2 = 0.0f;             // chi**2 of the fit
    m_quality = 0.0f;          // quality measure of the track

    m_firstLayer = 0;         // layer where the particle starts

    m_energyInput = 0.0f;         // energy input from the cal
    m_energyDetermined = 0.0f;    // energy determined using the tracker + cal
    
    m_hits = 0; 
    m_locator = 0; 

    Create();
}

TkrTrack::~TkrTrack() {
    Clean();
}

void TkrTrack::Clean() {

    if (m_hits) {
        int nEntries = m_hits->GetEntries();
        for (int i=0; i<nEntries; i++)
            delete m_hits->At(i);
        m_hits->Clear();
        delete m_hits;
        m_hits = 0;
    }

    if (m_locator) {
        int nEntries = m_locator->GetEntries();
        for (int i=0; i<nEntries; i++)
            delete m_locator->At(i);
        m_locator->Clear();
        delete m_locator;
        m_locator = 0;
    }

    m_id = 0;              // track id
    m_type = NONE;            // type of the track (gamma, e+, e-, hadron ...)
    m_motherId = NONE;        // id of the mother particle (-1 if none)

    m_residual = 0.0f;         // residual  
    m_chi2 = 0.0f;             // chi**2 of the fit
    m_quality = 0.0f;          // quality measure of the track

    m_firstLayer = 0;         // layer where the particle starts

    m_energyInput = 0.0f;         // energy input from the cal
    m_energyDetermined = 0.0f;    // energy determined using the tracker + cal
}

void TkrTrack::Create() {
    if (!m_hits) m_hits = new TObjArray();
    if (!m_locator) m_locator = new TObjArray();
}
