#include "reconRootData/TkrVecPointInfo.h"
#include "Riostream.h"

#include "ReconObjectManager.h"

ClassImp(TkrVecPointInfo)

void TkrVecPointInfo::Clear(Option_t* /* option */) {
    // Make sure everything is properly zeroed
    m_maxNumSkippedLayers    = 0;
    m_numVecPoints           = 0;
    m_numBiLayersWVecPoints  = 0;
    m_maxNumLinkCombinations = 0.;
}

void TkrVecPointInfo::Print(Option_t *option) const {
    TObject::Print(option);
}

//======================================================
// For Unit Tests
//======================================================

void TkrVecPointInfo::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrVecPointInfo::CompareInRange( const TkrVecPointInfo & /*ref*/, const std::string & name ) const {

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
