#include "reconRootData/TkrRecon.h"
#include "reconRootData/TkrTrack.h"
#include "reconRootData/TkrSiCluster.h"

ClassImp(TkrRecon)

TkrRecon::TkrRecon() {
    // default constructor
    m_tracks = 0;
//    m_gamma = 0;
    m_SiClusters = 0;
}


TkrRecon::~TkrRecon() {
    // destructor
    Clean();
}

void TkrRecon::Clean() {
    // Clears and destroys all TObjArrays contained in the TkrRecon object.
    // Necessary both for destruction and between subsequent fills and
    // writes/reads to/from file.
   
    if (m_tracks) {
        int nEntries = m_tracks->GetEntries();
        for (int i=0; i<nEntries; i++)
            delete m_tracks->At(i);
        m_tracks->Clear();
        delete m_tracks;
        m_tracks = 0;
    }

    //if (m_gamma) {
    //    delete m_gamma;
    //    m_gamma = 0;
   // }

    if (m_SiClusters) {
        int nEntries = m_SiClusters->GetEntries();
        for (int i=0; i<nEntries; i++)
            delete m_SiClusters->At(i);
        m_SiClusters->Clear();
        delete m_SiClusters;
        m_SiClusters = 0;
    }
}

void TkrRecon::Create() {
    if (!m_tracks) m_tracks = new TObjArray();
    if (!m_SiClusters) m_SiClusters = new TObjArray();
    m_gamma.Create();
}
