#include "reconRootData/TkrSiCluster.h"


ClassImp(TkrSiCluster)

TkrSiCluster::TkrSiCluster() {
    // default constructor
  m_id = 0;               // id of the cluster
  m_layer = 0;           // layer
  m_xy = X;             // x-y layer
  
  m_centerStrip = 0;    // center strip of the cluster  HMA ?:  why is this a float?
  m_numStrips = 0;     // number of strips    
  
  m_position = 0.0f;       // strip position  
  m_z = 0.0f;              // z position   (HMA Question: z position of the layer???)
}

TkrSiCluster::TkrSiCluster(UInt_t id) {
    // non-default constructor
  m_id = id;               // id of the cluster
  m_layer = 0;           // layer
  m_xy = X;             // x-y layer
  
  m_centerStrip = 0;    // center strip of the cluster  HMA ?:  why is this a float?
  m_numStrips = 0;     // number of strips    
  
  m_position = 0.0f;       // strip position  
  m_z = 0.0f;              // z position   
}


TkrSiCluster::~TkrSiCluster() {
    Clean();
}

void TkrSiCluster::Clean() {
  m_id = 0;               // id of the cluster
  m_layer = 0;           // layer
  m_xy = X;             // x-y layer
  
  m_centerStrip = 0;    // center strip of the cluster  HMA ?:  why is this a float?
  m_numStrips = 0;     // number of strips    
  
  m_position = 0.0f;       // strip position  
  m_z = 0.0f;              // z position   (HMA Question: z position of the layer???)
}
