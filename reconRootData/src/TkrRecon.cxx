#include "reconRootData/TkrRecon.h"
#include <iostream>

ClassImp(TkrRecon)

TkrRecon::TkrRecon() {
    m_clusterCol = 0;
    m_trackCol = 0;
    m_vertexCol = 0;
}

TkrRecon::~TkrRecon() {
    Clear();
    if (m_clusterCol) delete m_clusterCol;
    m_clusterCol = 0;
    if (m_trackCol) delete m_trackCol;
    m_trackCol = 0;
    if (m_vertexCol) delete m_vertexCol;
    m_vertexCol = 0;
}

void TkrRecon::initialize() {
    if (!m_clusterCol) m_clusterCol = new TObjArray();
    if (!m_trackCol) m_trackCol = new TObjArray();
    if (!m_vertexCol) m_vertexCol = new TObjArray();
}

void TkrRecon::Clear(Option_t* /* option */) {
    if (m_clusterCol) m_clusterCol->Delete();
    if (m_trackCol) m_trackCol->Delete();
    if (m_vertexCol) m_vertexCol->Delete();
}

void TkrRecon::Print(Option_t *option) const {
    TObject::Print(option);
    using namespace std;
    cout << "Number of TkrClusters: " << m_clusterCol->GetEntries() << endl;
    cout << "Number of Tracks: " << m_trackCol->GetEntries() << endl;
    cout << "Number of Vertices: " << m_vertexCol->GetEntries() << endl;
}
