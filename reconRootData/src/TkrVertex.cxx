#include "reconRootData/TkrVertex.h"
#include <iostream>

ClassImp(TkrVertex)

TkrVertex::TkrVertex() 
{
    Clear();
}
    
TkrVertex::~TkrVertex() {
    Clear();
}

void TkrVertex::initializeInfo(UInt_t id, UInt_t layer, UInt_t tower, 
                               Double_t quality, Double_t energy)
{
    m_id         = id;
    m_energy     = energy;
    m_quality    = quality;
    m_itower     = tower;
    m_firstLayer = layer;
}
    
void TkrVertex::initializeVals(const TkrParams& vtxPar, 
                               const TkrCovMat& vtxCov, const TVector3& pos, const TVector3& dir)
{
    m_vertexPar = vtxPar;
    m_vertexCov = vtxCov;
    m_position  = pos;
    m_direction = dir;
}

void TkrVertex::Clear(Option_t *option) {
    m_id          = -9999;
    m_energy      = 0.;
    m_quality     = -9999.;
    m_firstLayer  = -1;
    m_itower      = -1;

    m_vertexPar.Clear();
    //m_vertexCov.Clear();
    m_position = TVector3(0.,0.,0.);
    m_direction = TVector3(0., 0., 0.);

    m_trackIds.clear();
}

void TkrVertex::Print(Option_t *option) const {
    TObject::Print(option);
    using namespace std;
    cout << "Id: " << m_id << endl;
    cout << "Pos: (" << m_position.X() << "," << m_position.Y()
        << "," << m_position.Z() << ")  Dir: (" 
        << m_direction.X() << "," << m_direction.Y() << ","
        << m_direction.Z() << ")" << endl;
}
