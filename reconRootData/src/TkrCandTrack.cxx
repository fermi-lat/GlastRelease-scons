#include "reconRootData/TkrCandTrack.h"

ClassImp(TkrCandTrack)

TkrCandTrack::TkrCandTrack() {
    Clear();
}

TkrCandTrack::TkrCandTrack(UInt_t id, UInt_t layer, UInt_t tower, Double_t q, Double_t e, const TVector3& pos, const TVector3& dir)
{
    initialize(id, layer, tower, q, e, pos, dir);
}
    
void TkrCandTrack::initialize(UInt_t id, UInt_t layer, UInt_t tower, Double_t q, Double_t e, const TVector3& pos, const TVector3& dir)
{
    m_id          = id;
    m_energy      = e;
    m_quality     = q;
    m_firstLayer  = layer;
    m_itower      = tower;

    m_position    = pos;
    m_direction   = dir;

    m_hits.clear();
}
   
void TkrCandTrack::Clear(Option_t *option) {
    // Make sure everything is properly zeroed
    m_id          = -9999;
    m_energy      = 0.;
    m_quality     = -9999.;
    m_firstLayer  = -1;
    m_itower      = -1;

    m_position = TVector3(0.,0.,0.);
    m_direction = TVector3(0., 0., 0.);

    m_hits.clear();
}

void TkrCandTrack::Print(Option_t *option) const {
    TObject::Print(option);
}
