#include "reconRootData/TkrVertex.h"

ClassImp(TkrVertex)

// Default constructor
TkrVertex::TkrVertex() 
{
    // Make sure everything is properly zeroed
    m_id          = -9999.;
    m_energy      = 0.;
    m_quality     = -9999.;
    m_firstLayer  = -1;
    m_itower      = -1;

    m_vertexPar.Clear();
    m_vertexCov.Clear();
    m_position.Clear();
    m_direction.Clear();

    m_trackIds.clear();
}
    
void TkrVertex::initializeInfo(UInt_t id, UInt_t layer, UInt_t tower, Double_t quality, Double_t energy)
{
    m_id         = id;
    m_energy     = energy;
    m_quality    = quality;
    m_itower     = tower;
    m_firstLayer = layer;
}
    
void TkrVertex::initializeVals(const TkrParams& vtxPar, const TkrCovMat& vtxCov, const TVectorD& pos, const TVectorD& dir)
{
    m_vertexPar = vtxPar;
    m_vertexCov = vtxCov;
    m_position  = pos;
    m_direction = dir;
}
