#include "reconRootData/TkrFilterParams.h"
#include "Riostream.h"

#include "ReconObjectManager.h"

ClassImp(TkrFilterParams)
void TkrFilterParams::initializeInfo(TVector3& position, 
                                     TVector3& direction,
                                     UInt_t    statusBits,
                                     Double_t  energy,
                                     Int_t     numBiLayers,
                                     Int_t     numIterations,
                                     Int_t     numHitsTotal,
                                     Int_t     numDropped,
                                     Double_t  chiSquare,
                                     Double_t  aveDist,
                                     Double_t  transRms, 
                                     Double_t  longRms,
                                     Double_t  longRmsAsym
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
    m_aveDist        = aveDist;
    m_transRms       = transRms;
    m_longRms        = longRms;
    m_longRmsAsym    = longRmsAsym;
}

void TkrFilterParams::Clear(Option_t* /* option */) {
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
    m_aveDist        = 0.;
    m_transRms       = 0.;
    m_longRms        = 0.;
    m_longRmsAsym    = 0.;
}

void TkrFilterParams::Print(Option_t *option) const {
    TObject::Print(option);
}



//======================================================
// For Unit Tests
//======================================================

void TkrFilterParams::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrFilterParams::CompareInRange( const TkrFilterParams & /*ref*/, const std::string & name ) const {

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

void* TkrFilterParams::operator new(size_t /*size*/)
{
    TkrFilterParams* temp = ReconObjectManager::getPointer()->getNewTkrFilterParams();

    return temp;
}

void* TkrFilterParams::operator new(size_t /*size*/, void* vp)
{
    return vp;
}

void TkrFilterParams::operator delete(void* /*p*/)
{
    // Since we let ReconObjectManager handles memory, nothing to do here
    return;
}

