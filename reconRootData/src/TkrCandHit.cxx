#include "reconRootData/TkrCandHit.h"

ClassImp(TkrCandHit)

TkrCandHit::TkrCandHit()
{
    Clear();
}

TkrCandHit::TkrCandhit(const TVector3& pos, UInt_t hitId, UInt_t tower, UInt_t layer, AXIS view)
{
    initialize(pos,hitId,tower,layer,view);
}

void TkrCandHit::initialize(const TVector3& pos, UInt_t hitId, UInt_t tower, UInt_t layer, AXIS view)
{
    m_position = pos;
    m_hitIndex = hitId;
    m_tower    = tower;
    m_plane    = layer;
    m_view     = view;
}

void TkrCandHit::Clear(Option_t *option) {
    m_position = TVector3(0., 0., 0.);
    m_hitIndex = -1;
    m_tower    = -1;
    m_plane    = -1;
    m_view     = NONE;
}

void TkrCandHit::Print(Option_t *option) const {
    TObject::Print(option);
}