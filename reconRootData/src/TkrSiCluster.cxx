#include "reconRootData/TkrSiCluster.h"


ClassImp(TkrSiCluster)

TkrSiCluster::TkrSiCluster() {
    Clear();
}

TkrSiCluster::TkrSiCluster(UInt_t id) {
    Clear();
    m_id = id;
}


TkrSiCluster::~TkrSiCluster() {
    Clear();
}

void TkrSiCluster::Clear() {
    m_id = 0;               
    m_layer = 0;
    m_xy = X;
    
    m_centerStrip = 0; 
    m_numStrips = 0;      
    
    m_position = 0.0f; 
    m_z = 0.0f;
}
