
#include "digiRootData/TowerId.h"

ClassImp(TowerId)

TowerId::TowerId (UShort_t id):m_id(id)
{
}

TowerId::TowerId (UInt_t ix, UInt_t iy):m_id(ix+xNum*iy) 
{
}


Bool_t TowerId::neighbor (const TowerId& n) const
{
    UInt_t dx = ix()-n.ix(), dy=iy()-n.iy();
    return dx<2 && dx>-2 && dy<2 && dy>-2; 
}

