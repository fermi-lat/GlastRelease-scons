
#include "reconRootData/AdfRecon.h"
#include "commonRootData/RootDataUtil.h"
#include "Riostream.h"
#include "TCollection.h"  // Declares TIter

using namespace std;

ClassImp(reconRootData::AdfRecon)

namespace reconRootData {

AdfRecon::AdfRecon() {
    m_taggerClusterCol = 0; 
    m_qdcHitCol = 0; 
    m_numQdcHit = -1;

    m_eventNumber = 0;
    m_spillNumber = 0;
}


AdfRecon::~AdfRecon() {
  
    if (m_taggerClusterCol) {
        m_taggerClusterCol->Delete();
        delete m_taggerClusterCol;
        m_taggerClusterCol = 0;
    }
  
    if (m_qdcHitCol) {
        m_qdcHitCol->Delete();
        delete m_qdcHitCol;
        m_qdcHitCol = 0;
    }
    
}



void AdfRecon::Clear(Option_t *option) {

    m_eventNumber = 0;
    m_spillNumber = 0;

    if (m_taggerClusterCol) m_taggerClusterCol->Clear("C");
    if (m_qdcHitCol) m_qdcHitCol->Clear("C");
    m_numQdcHit = -1;

}

void AdfRecon::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "AdfRecon:" << endl;
    cout << "EventNumber: " << m_eventNumber << " SpillNumber: " 
         << m_spillNumber << endl;
    cout << "Num Tagger Clusters: " << m_taggerClusterCol->GetEntriesFast() 
         << " Num Qdc Hits: " << m_numQdcHit << endl;
    cout << dec;
}

void  AdfRecon::addTaggerCluster(TaggerCluster* cluster) { 
    if (!m_taggerClusterCol) m_taggerClusterCol = new TObjArray();
    m_taggerClusterCol->Add(cluster);
}

const reconRootData::TaggerCluster* AdfRecon::getTaggerCluster(UInt_t ind) const {
    if (((Int_t)ind) < m_taggerClusterCol->GetEntriesFast()) 
        return ((reconRootData::TaggerCluster*)m_taggerClusterCol->At(ind));
    else
        return 0;
}

commonRootData::QdcHit* AdfRecon::addQdcHit(UInt_t channel, UInt_t pulseHgt, Bool_t isPedSub) {
    // Add a new QdcHit entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    if (!m_qdcHitCol) m_qdcHitCol = new TClonesArray("commonRootData::QdcHit",1);
    ++m_numQdcHit;
    TClonesArray &localCol = *m_qdcHitCol;
    new(localCol[m_numQdcHit]) commonRootData::QdcHit(channel, pulseHgt, isPedSub);
    return ((commonRootData::QdcHit*)(localCol[m_numQdcHit]));
}

const commonRootData::QdcHit* AdfRecon::getQdcHit(UInt_t ind) const {
    if (((Int_t)ind) < m_qdcHitCol->GetEntriesFast()) 
        return ((commonRootData::QdcHit*)m_qdcHitCol->At(ind));
    else
        return 0;
}

void AdfRecon::Fake( Int_t ievent, Float_t randNum ) {
    setEventNumber(ievent);
    setSpillNumber(5);
    // Add 2 TaggerClusters
    reconRootData::TaggerCluster* clus1 = new TaggerCluster();
    reconRootData::TaggerCluster* clus2 = new TaggerCluster();

    clus1->Fake(ievent, randNum);
    clus2->Fake(ievent, randNum);

    addTaggerCluster(clus1);
    addTaggerCluster(clus2);
   
    // Add 3 QdcHits
    addQdcHit(1,2,true);
    addQdcHit(2,3,false);
    addQdcHit(3,4,true);
}

Bool_t AdfRecon::CompareInRange( const AdfRecon &ref, const std::string& name ) const {
    bool result = true;

    result = rootdatautil::CompareInRange(getEventNumber(),ref.getEventNumber(),"EventNumber") && result;
    result = rootdatautil::CompareInRange(getSpillNumber(),ref.getSpillNumber(),"SpillNumber") && result;

    result = rootdatautil::TObjArrayCompareInRange<reconRootData::TaggerCluster>(m_taggerClusterCol,ref.m_taggerClusterCol) && result;
    result = rootdatautil::TObjArrayCompareInRange<commonRootData::QdcHit>(m_qdcHitCol,ref.m_qdcHitCol) && result;

    if (!result) {
        if (name == "") {
            std::cout<<"Comparison ERROR for " << ClassName() << std::endl;
        } else {
            std::cout << "Comparison ERROR for " << name<<std::endl;
        }
    }

    return result ;

}

} // end namespace
