#include "reconRootData/TkrRecon.h"
#include <iostream>

ClassImp(TkrRecon)

TkrRecon::TkrRecon() {
    m_siClusterCol = 0;
    m_trackCandCol = 0;
    m_trackCol = 0;
    m_vertexCol = 0;
}

TkrRecon::~TkrRecon() {
    Clear();
    if (m_siClusterCol) delete m_siClusterCol;
    m_siClusterCol = 0;
    if (m_trackCandCol) delete m_trackCandCol;
    m_trackCandCol = 0;
    if (m_trackCol) delete m_trackCol;
    m_trackCol = 0;
    if (m_vertexCol) delete m_vertexCol;
    m_vertexCol = 0;
}

void TkrRecon::initialize() {
    if (!m_siClusterCol) m_siClusterCol = new TObjArray();
    if (!m_trackCandCol) m_trackCandCol = new TObjArray();
    if (!m_trackCol) m_trackCol = new TObjArray();
    if (!m_vertexCol) m_vertexCol = new TObjArray();
}

void TkrRecon::Clear(Option_t *option) {
    if (m_siClusterCol) m_siClusterCol->Delete();
    if (m_trackCandCol) m_trackCandCol->Delete();
    if (m_trackCol) m_trackCol->Delete();
    if (m_vertexCol) m_vertexCol->Delete();
}

void TkrRecon::Print(Option_t *option) const {
    TObject::Print(option);
    using namespace std;
    cout << "Number of TkrSiClusters: " << m_siClusterCol->GetEntries() << endl;
    cout << "Number of Track Candidates: " << m_trackCandCol->GetEntries() << endl;
    cout << "Number of Tracks: " << m_trackCol->GetEntries() << endl;
    cout << "Number of Vertics: " << m_vertexCol->GetEntries() << endl;
}