#include "commonRootData/idents/TowerId.h"
#include <iostream>

ClassImp(commonRootData::TowerId)

namespace commonRootData {

TowerId::TowerId (UShort_t id):m_id(id)
{
}

TowerId::TowerId (UInt_t ix, UInt_t iy):m_id(ix+xNum*iy) 
{
}

void TowerId::Clear(Option_t *option) {

}

void TowerId::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Id: " << m_id << " (x, y): (" << ix() << ", " << iy() << ")" << endl;
}


Bool_t TowerId::neighbor (const TowerId& n) const
{
    Int_t dx = ix()-n.ix(), dy=iy()-n.iy();
    return dx<2 && dx>-2 && dy<2 && dy>-2; 
}

} // end namespace

