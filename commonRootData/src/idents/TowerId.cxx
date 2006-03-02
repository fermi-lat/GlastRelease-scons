#include "commonRootData/idents/TowerId.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(TowerId)

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


#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)
#define COMPARE_IN_RANGE_NO_GET(att) rootdatautil::CompareInRange(att(),ref. ## att(),#att)

Bool_t TowerId::CompareInRange( const TowerId & ref, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_IN_RANGE_NO_GET(id) && result ;

    if (!result) {
        std::cout<<"Comparison ERROR for "<<name<<std::endl ;
    }
    return result ;

}

