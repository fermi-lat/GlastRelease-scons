#include "reconRootData/TkrSiCluster.h"
#include <iostream>


ClassImp(TkrSiCluster)

TkrSiCluster::TkrSiCluster() {
    Clear();
}

TkrSiCluster::~TkrSiCluster() {
    Clear();
}


void TkrSiCluster::initialize(UInt_t id, UInt_t layer, TKRAxes xy, UShort_t center,
                         UShort_t numStrips, Float_t pos, Float_t zpos) {
    Clear();
    m_id = id;
    m_layer = layer;
    m_xy = xy;
    m_centerStrip = center;
    m_numStrips = numStrips;
    m_position = pos;
    m_z = zpos;
}

void TkrSiCluster::Clear(Option_t *option) {
    m_id = 0;               
    m_layer = 0;
    m_xy = X;
    
    m_centerStrip = 0; 
    m_numStrips = 0;      
    
    m_position = 0.0f; 
    m_z = 0.0f;
}

void TkrSiCluster::Print(Option_t *option) const {
    TObject::Print(option);
    using namespace std;
    cout << "Id: " << m_id << " Layer: " << m_layer 
        << " Axes: " << m_xy << endl;
    cout << "Center, NumStrips: " << m_centerStrip << ","
        << m_numStrips << endl;
    cout << "Pos, ZPos: " << m_position << "," << m_z << endl;
}
