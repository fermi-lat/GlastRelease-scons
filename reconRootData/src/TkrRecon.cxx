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
    const Int_t nd = 20000;
    static Int_t indCluster = 0, indTrack=0, indVertex=0;
    static Int_t limitCluster = 100;
    static Int_t limitTrack = 100;
    static Int_t limitVertex = 100;
    static TkrCluster *keepCluster[nd];
    static TkrVertex *keepVertex[nd];
    static TkrTrack *keepTrack[nd];

//    if (m_clusterCol) m_clusterCol->Delete();
//    if (m_trackCol) m_trackCol->Delete();
//    if (m_vertexCol) m_vertexCol->Delete();

    // Emulate what the TClonesArray is doing as introduced in DigiEvent
    // For the Cluster TObjArray, by deleting the objects in chunks
    Int_t i, j;
    Int_t n = m_clusterCol->GetEntries();
    if (n > limitCluster) {
        limitCluster = n + 10;
        if (limitCluster > nd)  
            std::cout << "!!!Warning: limit for TkrCluster Array is greater than " 
                 << nd << std::endl;
        for (j = 0; j<indCluster; j++) delete keepCluster[j];
        indCluster = 0;
    }
    for (i=0;i<n;i++) 
        keepCluster[indCluster+i] = (TkrCluster*)m_clusterCol->At(i);
    indCluster +=n;
    if (indCluster > nd-limitCluster) {
        for (j = 0; j<indCluster;j++) delete keepCluster[j];
        indCluster = 0;
    }
    m_clusterCol->Clear();

    // Emulate what the TClonesArray is doing as introduced in DigiEvent
    // For the Track TObjArray, by deleting the objects in chunks
    n = m_trackCol->GetEntries();
    if (n > limitTrack) {
        limitTrack = n + 10;
        if (limitTrack > nd)  
            std::cout << "!!!Warning: limit for TkrTrack Array is greater than " 
                 << nd << std::endl;
        for (j = 0; j<indTrack; j++) delete keepTrack[j];
        indTrack = 0;
    }
    for (i=0;i<n;i++) keepTrack[indTrack+i] = (TkrTrack*)m_trackCol->At(i);
    indTrack +=n;
    if (indTrack > nd-limitTrack) {
        for (j = 0; j<indTrack;j++) delete keepTrack[j];
        indTrack = 0;
    }
    m_trackCol->Clear();

    // Emulate what the TClonesArray is doing as introduced in DigiEvent
    n = m_vertexCol->GetEntries();
    if (n > limitVertex) {
        limitVertex = n + 10;
        if (limitVertex > nd)  
            std::cout << "!!!Warning: limit for TkrVertex Array is greater than " 
                 << nd << std::endl;
        for (j = 0; j<indVertex; j++) delete keepVertex[j];
        indVertex = 0;
    }
    for (i=0;i<n;i++) keepVertex[indVertex+i] = (TkrVertex*)m_vertexCol->At(i);
    indVertex +=n;
    if (indVertex > nd-limitVertex) {
        for (j = 0; j<indVertex;j++) delete keepVertex[j];
        indVertex = 0;
    }
    m_vertexCol->Clear();
}

void TkrRecon::Print(Option_t *option) const {
    TObject::Print(option);
    using namespace std;
    cout << "Number of TkrClusters: " << m_clusterCol->GetEntries() << endl;
    cout << "Number of Tracks: " << m_trackCol->GetEntries() << endl;
    cout << "Number of Vertices: " << m_vertexCol->GetEntries() << endl;
}
