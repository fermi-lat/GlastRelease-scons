#include "reconRootData/TkrCandTrack.h"

ClassImp(TkrCandTrack)

// Default constructor
TkrCandTrack::TkrCandTrack() 
{
    // Make sure everything is properly zeroed
    m_id          = -9999;
    m_energy      = 0.;
    m_quality     = -9999.;
    m_firstLayer  = -1;
    m_itower      = -1;

    m_position.Zero();
    m_direction.Zero();

    m_hits.clear();
}

TkrCandTrack::TkrCandTrack(UInt_t id, UInt_t layer, UInt_t tower, Double_t q, Double_t e, const TVectorD& pos, const TVectorD& dir)
{
    initialize(id, layer, tower, q, e, pos, dir);
}
    
void TkrCandTrack::initialize(UInt_t id, UInt_t layer, UInt_t tower, Double_t q, Double_t e, const TVectorD& pos, const TVectorD& dir)
{
    m_id          = id;
    m_energy      = e;
    m_quality     = q;
    m_firstLayer  = layer;
    m_itower      = tower;

    m_position.ResizeTo(pos);
    m_position    = pos;
    m_direction.ResizeTo(dir);
    m_direction   = dir;

    m_hits.clear();
}
   
