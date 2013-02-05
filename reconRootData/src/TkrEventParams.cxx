#include "reconRootData/TkrEventParams.h"
#include "Riostream.h"

#include "ReconObjectManager.h"

ClassImp(TkrEventParams)
void TkrEventParams::initializeInfo(TVector3& position, 
                                     TVector3& direction,
                                     UInt_t    statusBits,
                                     Double_t  energy,
                                     Int_t     numBiLayers,
                                     Int_t     numIterations,
                                     Int_t     numHitsTotal,
                                     Int_t     numDropped,
                                     Double_t  chiSquare,
                                     Double_t  transRms, 
                                     Double_t  longRms
                                     )
{
    m_statusBits     = statusBits;
    m_eventEnergy    = energy;
    m_eventPosition  = position;
    m_eventAxis      = direction;
    m_numBiLayers    = numBiLayers;
    m_numIterations  = numIterations;
    m_numHitsTotal   = numHitsTotal;
    m_numDropped     = numDropped;
    m_chiSquare      = chiSquare;
    m_transRms       = transRms;
    m_longRmsAve     = longRms;
}

void TkrEventParams::Clear(Option_t* /* option */) {
    // Make sure everything is properly zeroed
    m_statusBits     = 0;
    m_eventEnergy    = 0.;
    m_eventPosition  = TVector3(0.,0.,0.);
    m_eventAxis      = TVector3(0.,0.,0.);
    m_numBiLayers    = 0;
    m_numIterations  = 0;
    m_numHitsTotal   = 0;
    m_numDropped     = 0;
    m_chiSquare      = -1.;
    m_transRms       = 0.;
    m_longRmsAve     = 0.;
}

void TkrEventParams::Print(Option_t *option) const {
    TObject::Print(option);
}



//======================================================
// For Unit Tests
//======================================================

void TkrEventParams::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrEventParams::CompareInRange( const TkrEventParams & /*ref*/, const std::string & name ) const {

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
