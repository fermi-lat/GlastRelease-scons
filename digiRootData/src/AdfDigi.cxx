
#include "digiRootData/AdfDigi.h"
#include "commonRootData/RootDataUtil.h"
#include "Riostream.h"
#include "TCollection.h"  // Declares TIter

using namespace std;

ClassImp(AdfDigi)

AdfDigi::AdfDigi() {
    m_taggerHitCol = 0; 
    m_numTaggerHit = -1;
    m_qdcHitCol = 0; 
    m_numQdcHit = -1;
    m_scalerHitCol = 0;
    m_numScalerHit = 0;

    m_eventNumber = 0;
    m_spillNumber = 0;
}


AdfDigi::~AdfDigi() {
  
    if (m_taggerHitCol) {
        m_taggerHitCol->Delete();
        delete m_taggerHitCol;
        m_taggerHitCol = 0;
    }
  
    if (m_qdcHitCol) {
        m_qdcHitCol->Delete();
        delete m_qdcHitCol;
        m_qdcHitCol = 0;
    }
    

    if (m_scalerHitCol) {
        m_scalerHitCol->Delete();
        delete m_scalerHitCol;
        m_scalerHitCol = 0;
    }
}



void AdfDigi::Clear(Option_t *option) {

    TObject::Clear(option);
    m_eventNumber = 0;
    m_spillNumber = 0;

    if (m_taggerHitCol) m_taggerHitCol->Clear("C");
    if (m_qdcHitCol) m_qdcHitCol->Clear("C");
    if (m_scalerHitCol) m_scalerHitCol->Clear("C");
    m_numTaggerHit = -1;
    m_numQdcHit = -1;
    m_numScalerHit = -1;

}

void AdfDigi::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "AdfDigi:" << endl;
    cout << "EventNumber: " << m_eventNumber << " SpillNumber: " 
         << m_spillNumber << endl;
    cout << "Num Tagger Hits: " << m_numTaggerHit+1 << " Num Qdc Hits: " 
         << m_numQdcHit+1 << " Num Scaler Hits: " << m_numScalerHit+1 << endl;
    cout << dec;
}

commonRootData::TaggerHit* AdfDigi::addTaggerHit(UInt_t moduleId, UInt_t layerId, 
                                 UInt_t stripId, Double_t pulseHgt, 
                                 Double_t sigma, Bool_t isPedSubtracted) {
    // Add a new TaggerHit entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    if (!m_taggerHitCol) m_taggerHitCol = new TClonesArray("commonRootData::TaggerHit", 1);
    ++m_numTaggerHit;
    TClonesArray &localCol = *m_taggerHitCol;
    new(localCol[m_numTaggerHit]) commonRootData::TaggerHit(moduleId, layerId, stripId, pulseHgt, sigma, isPedSubtracted);
    return ((commonRootData::TaggerHit*)(localCol[m_numTaggerHit]));
}

const commonRootData::TaggerHit* AdfDigi::getTaggerHit(UInt_t ind) const {
    if (((Int_t)ind) < m_taggerHitCol->GetEntriesFast()) 
        return ((commonRootData::TaggerHit*)m_taggerHitCol->At(ind));
    else
        return 0;
}

commonRootData::QdcHit* AdfDigi::addQdcHit(UInt_t channel, Double_t pulseHgt, UInt_t mod, Double_t sigma, Bool_t isPedSub) {
    // Add a new QdcHit entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    if (!m_qdcHitCol) m_qdcHitCol = new TClonesArray("commonRootData::QdcHit",1);
    ++m_numQdcHit;
    TClonesArray &localCol = *m_qdcHitCol;
    new(localCol[m_numQdcHit]) commonRootData::QdcHit(channel, pulseHgt, mod, sigma, isPedSub);
    return ((commonRootData::QdcHit*)(localCol[m_numQdcHit]));
}

const commonRootData::QdcHit* AdfDigi::getQdcHit(UInt_t ind) const {
    if (((Int_t)ind) < m_qdcHitCol->GetEntriesFast()) 
        return ((commonRootData::QdcHit*)m_qdcHitCol->At(ind));
    else
        return 0;
}

commonRootData::ScalerHit* AdfDigi::addScalerHit(UInt_t channel, UInt_t val) {
    // Add a new ScalerHit entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    if (!m_scalerHitCol) m_scalerHitCol = new TClonesArray("commonRootData::ScalerHit",1);
    ++m_numScalerHit;
    TClonesArray &localCol = *m_scalerHitCol;
    new(localCol[m_numScalerHit]) commonRootData::ScalerHit(channel,val);
    return ((commonRootData::ScalerHit*)(localCol[m_numScalerHit]));
}

const commonRootData::ScalerHit* AdfDigi::getScalerHit(UInt_t ind) const {
    if (((Int_t)ind) < m_scalerHitCol->GetEntriesFast()) 
        return ((commonRootData::ScalerHit*)m_scalerHitCol->At(ind));
    else
        return 0;
}


void AdfDigi::Fake( Int_t ievent, Float_t /*randNum*/ ) {
    setEventNumber(ievent);
    setSpillNumber(5);
    // Add 2 TaggerHits
    addTaggerHit(1,2,3,4.4,5.5,true);
    addTaggerHit(5,6,7,8.8,9.9,false);
   
    // Add 3 QdcHits
    addQdcHit(1,2.2,3,4.4,true);
    addQdcHit(2,3.3,4,5.5,false);
    addQdcHit(3,4.4,5,6.6,true);

    addScalerHit(0,1);
    addScalerHit(1,2);
}

Bool_t AdfDigi::CompareInRange( const AdfDigi &ref, const std::string& name ) const {
    bool result = true;

    result = rootdatautil::CompareInRange(getEventNumber(),ref.getEventNumber(),"EventNumber") && result;
    result = rootdatautil::CompareInRange(getSpillNumber(),ref.getSpillNumber(),"SpillNumber") && result;

    result = rootdatautil::TObjArrayCompareInRange<commonRootData::TaggerHit>(m_taggerHitCol,ref.m_taggerHitCol) && result;
    result = rootdatautil::TObjArrayCompareInRange<commonRootData::QdcHit>(m_qdcHitCol,ref.m_qdcHitCol) && result;
    result = rootdatautil::TObjArrayCompareInRange<commonRootData::ScalerHit>(m_scalerHitCol,ref.m_scalerHitCol) && result;

    if (!result) {
        if (name == "") {
            std::cout<<"Comparison ERROR for " << ClassName() << std::endl;
        } else {
            std::cout << "Comparison ERROR for " << name<<std::endl;
        }
    }

    return result ;

}
