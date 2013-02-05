#include "reconRootData/TkrVecPoint.h"
#include "Riostream.h"

#include "ReconObjectManager.h"

ClassImp(TkrVecPoint)

void TkrVecPoint::Clear(Option_t* /* option */) {
    // Make sure everything is properly zeroed
    m_status    = 0;
    m_layer     = 0;
    m_pXCluster = 0;
    m_pYCluster = 0;
    m_position  = TVector3(0.,0.,0.);

}

void TkrVecPoint::Print(Option_t *option) const {
    TObject::Print(option);
}



//======================================================
// For Unit Tests
//======================================================

void TkrVecPoint::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrVecPoint::CompareInRange( const TkrVecPoint & /*ref*/, const std::string & name ) const {

    bool result = true ;

    //result = COMPARE_IN_RANGE(Energy) && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}

void* TkrVecPoint::operator new(size_t /*size*/)
{
    TkrVecPoint* temp = ReconObjectManager::getPointer()->getNewTkrVecPoint();

    return temp;
}

void* TkrVecPoint::operator new(size_t /*size*/, void* vp)
{
    return vp;
}

void TkrVecPoint::operator delete(void* /*p*/)
{
    // Since we let ReconObjectManager handles memory, nothing to do here
    return;
}

