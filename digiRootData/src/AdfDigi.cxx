//
// The Event class is the top level class for the ASCII to ROOT
// conversion. It contains all the information about a single event.
//

#include "digiRootData/AdfDigi.h"
#include "commonRootData/RootDataUtil.h"
#include "Riostream.h"
#include "TCollection.h"  // Declares TIter

using namespace std;

ClassImp(AdfDigi)

// Allocate the TClonesArrays just once
TClonesArray *AdfDigi::s_taggerHitStaticCol = 0;
TClonesArray *AdfDigi::s_qdcHitStaticCol = 0;

AdfDigi::AdfDigi() {
    if (!s_taggerHitStaticCol) 
        s_taggerHitStaticCol = new TClonesArray("TaggerHit",1);
    m_taggerHitCol = s_taggerHitStaticCol;
    m_numTaggerHit = -1;

    // Assign default values to members
    if (!s_qdcHitStaticCol) s_qdcHitStaticCol = new TClonesArray("QdcHit", 1);
    m_qdcHitCol = s_qdcHitStaticCol;
    m_numQdcHit = -1;

    Clear();
}

AdfDigi::AdfDigi(const AdfDigi& copy):TObject(copy) {

    m_eventNumber = copy.m_eventNumber;
    m_spillNumber = copy.m_spillNumber;
    m_numTaggerHit = copy.m_numTaggerHit;
    m_numQdcHit = copy.m_numQdcHit;
    m_taggerHitCol = copy.m_taggerHitCol;
    m_qdcHitCol = copy.m_qdcHitCol;
}


AdfDigi::~AdfDigi() {
  
    if(m_taggerHitCol == s_taggerHitStaticCol) s_taggerHitStaticCol = 0;
    m_taggerHitCol->Delete();
    delete m_taggerHitCol;
    m_taggerHitCol = 0;
  
    if(m_qdcHitCol == s_qdcHitStaticCol) s_qdcHitStaticCol = 0;
    m_qdcHitCol->Delete();
    delete m_qdcHitCol;
    m_qdcHitCol = 0;
    
}

AdfDigi& AdfDigi::operator=(const AdfDigi& copy) {

    m_eventNumber = copy.m_eventNumber;
    m_spillNumber = copy.m_spillNumber;
    m_numTaggerHit = copy.m_numTaggerHit;
    m_numQdcHit = copy.m_numQdcHit;
    m_taggerHitCol = copy.m_taggerHitCol;
    m_qdcHitCol = copy.m_qdcHitCol;
    return *this;
}


void AdfDigi::Clear(Option_t *option) {

    m_eventNumber = 0;
    m_spillNumber = 0;

    m_taggerHitCol->Clear("C");
    m_qdcHitCol->Clear("C");
    m_numTaggerHit = -1;
    m_numQdcHit = -1;

}

void AdfDigi::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "AdfDigi:" << endl;
    cout << "EventNumber: " << m_eventNumber << " SpillNumber: " 
         << m_spillNumber << endl;
    cout << "Num Tagger Hits: " << m_numTaggerHit << " Num Qdc Hits: " 
         << m_numQdcHit << endl;
    cout << dec;
}

TaggerHit* AdfDigi::addTaggerHit(UInt_t moduleId, UInt_t layerId, 
                                 UInt_t stripId, UInt_t pulseHgt, 
                                 Bool_t isPedSubtracted) {
    // Add a new TaggerHit entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numTaggerHit;
    TClonesArray &localCol = *m_taggerHitCol;
    new(localCol[m_numTaggerHit]) TaggerHit(moduleId, layerId, stripId, 
                                            pulseHgt, isPedSubtracted);
    return ((TaggerHit*)(localCol[m_numTaggerHit]));
}

const TaggerHit* AdfDigi::getTaggerHit(UInt_t ind) const {
    if (((Int_t)ind) < m_taggerHitCol->GetEntriesFast()) 
        return ((TaggerHit*)m_taggerHitCol->At(ind));
    else
        return 0;
}

QdcHit* AdfDigi::addQdcHit(UInt_t channel, UInt_t pulseHgt, Bool_t isPedSub) {
    // Add a new QdcHit entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numQdcHit;
    TClonesArray &localCol = *m_qdcHitCol;
    new(localCol[m_numQdcHit]) QdcHit(channel, pulseHgt, isPedSub);
    return ((QdcHit*)(localCol[m_numQdcHit]));
}

const QdcHit* AdfDigi::getQdcHit(UInt_t ind) const {
    if (((Int_t)ind) < m_qdcHitCol->GetEntriesFast()) 
        return ((QdcHit*)m_qdcHitCol->At(ind));
    else
        return 0;
}

void AdfDigi::Fake( Int_t ievent, Float_t randNum ) {
    setEventNumber(ievent);
    setSpillNumber(5);
    // Add 2 TaggerHits
    addTaggerHit(1,2,3,4,true);
    addTaggerHit(5,6,7,8,false);
   
    // Add 3 QdcHits
    addQdcHit(1,2,true);
    addQdcHit(2,3,false);
    addQdcHit(3,4,true);
}

Bool_t AdfDigi::CompareInRange( AdfDigi &ref, const std::string& name ) {
    bool result = true;

    result = rootdatautil::CompareInRange(getEventNumber(),ref.getEventNumber(),"EventNumber") && result;
    result = rootdatautil::CompareInRange(getSpillNumber(),ref.getSpillNumber(),"SpillNumber") && result;

    result = rootdatautil::TObjArrayCompareInRange<TaggerHit>(m_taggerHitCol,ref.m_taggerHitCol) && result;
    result = rootdatautil::TObjArrayCompareInRange<QdcHit>(m_qdcHitCol,ref.m_qdcHitCol) && result;

    if (!result) {
        if (name == "") {
            std::cout<<"Comparison ERROR for " << ClassName() << std::endl;
        } else {
            std::cout << "Comparison ERROR for " << name<<std::endl;
        }
    }

    return result ;

}
