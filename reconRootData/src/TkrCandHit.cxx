#include "reconRootData/TkrCandHit.h"

ClassImp(TkrCandHit)

TkrCandHit::TkrCandHit()
{
    m_position.Zero();
    m_hitIndex = -1;
    m_tower    = -1;
    m_plane    = -1;
    m_view     = NONE;
}

TkrCandHit::TkrCandhit(const TVectorD& pos, UInt_t hitId, UInt_t tower, UInt_t layer, AXIS view)
{
    initialize(pos,hitId,tower,layer,view);
}

void TkrCandHit::initialize(const TVectorD& pos, UInt_t hitId, UInt_t tower, UInt_t layer, AXIS view)
{
    m_position.ResizeTo(pos);
    m_position = pos;
    m_hitIndex = hitId;
    m_tower    = tower;
    m_plane    = layer;
    m_view     = view;
}
