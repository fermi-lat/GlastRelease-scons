
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
    m_scalerHitCol = 0;
    m_numScalerHit = -1;

    m_eventNumber = 0;
    m_spillNumber = 0;
    m_px = 0.0; m_py = 0.0; m_pz = 0.0;
    m_eRec = 0.0; m_eCorr = 0.0;
    m_phiIn = 0.0; m_phiOut = 0.0;
    m_theta = 0.0; m_dphi = 0.0;
    m_numHighestCluster  = 0;
    m_nModule = 0;
    m_x = 0;
    m_y = 0;
    m_z = 0;
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
    
    if (m_scalerHitCol) {
        m_scalerHitCol->Delete();
        delete m_scalerHitCol;
        m_scalerHitCol = 0;
    }
    /*if (m_x) { delete [] m_x; m_x = 0; }
    if (m_y) { delete [] m_y; m_y = 0; }
    if (m_z) { delete [] m_z; m_z = 0; }*/
}



void AdfRecon::Clear(Option_t *option) {

    m_eventNumber = 0;
    m_spillNumber = 0;

    if (m_taggerClusterCol) m_taggerClusterCol->Clear("C");
    if (m_qdcHitCol) m_qdcHitCol->Clear("C");
    m_numQdcHit = -1;
    if (m_scalerHitCol) m_scalerHitCol->Clear("C");
    m_numScalerHit = -1;

    if (m_x) { delete [] m_x; m_x = 0; }
    if (m_y) { delete [] m_y; m_y = 0; }
    if (m_z) { delete [] m_z; m_z = 0;} 

    m_px = 0.0; m_py = 0.0; m_pz = 0.0;
    m_eRec = 0.0; m_eCorr = 0.0;
    m_phiIn = 0.0; m_phiOut = 0.0;
    m_theta = 0.0; m_dphi = 0.0;
    m_numHighestCluster  = 0;
    m_nModule = 0;

}

void AdfRecon::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "AdfRecon:" << endl;
    cout << "EventNumber: " << m_eventNumber << " SpillNumber: " 
         << m_spillNumber << endl;
    cout << "Num Tagger Clusters: " << m_taggerClusterCol->GetEntriesFast() 
         << "Num Scaler Hits: " << m_numScalerHit+1
         << " Num Qdc Hits: " << m_numQdcHit+1 << endl;
    cout << dec;
}

void AdfRecon::initXyz(const Double_t *x, const Double_t *y, const Double_t *z, UInt_t n) {
    m_nModule = n;
    m_x = new Double32_t[n];
    m_y = new Double32_t[n];
    m_z = new Double32_t[n];
    UInt_t i;
    for (i=0; i<n; i++) {
        m_x[i] = x[i];
        m_y[i] = y[i];
        m_z[i] = z[i];
    }
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

commonRootData::QdcHit* AdfRecon::addQdcHit(UInt_t channel, Double_t pulseHgt, UInt_t module, Double_t sigma, Bool_t isPedSub) {
    // Add a new QdcHit entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    if (!m_qdcHitCol) m_qdcHitCol = new TClonesArray("commonRootData::QdcHit",1);
    ++m_numQdcHit;
    TClonesArray &localCol = *m_qdcHitCol;
    new(localCol[m_numQdcHit]) commonRootData::QdcHit(channel, pulseHgt, module, sigma, isPedSub);
    return ((commonRootData::QdcHit*)(localCol[m_numQdcHit]));
}

const commonRootData::QdcHit* AdfRecon::getQdcHit(UInt_t ind) const {
    if (((Int_t)ind) < m_qdcHitCol->GetEntriesFast()) 
        return ((commonRootData::QdcHit*)m_qdcHitCol->At(ind));
    else
        return 0;
}

commonRootData::ScalerHit* AdfRecon::addScalerHit(UInt_t channel, UInt_t val) {
    // Add a new ScalerHit entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    if (!m_scalerHitCol) m_scalerHitCol = new TClonesArray("commonRootData::ScalerHit",1);
    ++m_numScalerHit;
    TClonesArray &localCol = *m_scalerHitCol;
    new(localCol[m_numScalerHit]) commonRootData::ScalerHit(channel, val);
    return ((commonRootData::ScalerHit*)(localCol[m_numScalerHit]));
}

const commonRootData::ScalerHit* AdfRecon::getScalerHit(UInt_t ind) const {
    if (((Int_t)ind) < m_scalerHitCol->GetEntriesFast()) 
        return ((commonRootData::ScalerHit*)m_scalerHitCol->At(ind));
    else
        return 0;
}
void AdfRecon::Fake( Int_t ievent, Float_t randNum ) {
    initEventNumber(ievent);
    initSpillNumber(5);
    // Add 2 TaggerClusters
    reconRootData::TaggerCluster* clus1 = new TaggerCluster();
    reconRootData::TaggerCluster* clus2 = new TaggerCluster();

    clus1->Fake(ievent, randNum);
    clus2->Fake(ievent, randNum);

    addTaggerCluster(clus1);
    addTaggerCluster(clus2);
   
    // Add 3 QdcHits
    addQdcHit(1,2.2,3,4.4,true);
    addQdcHit(2,3.3,4,5.5,false);
    addQdcHit(3,4.4,5,6.6,true);

    addScalerHit(6,7);

    initNumHighestCluster(8);
    init(1.1, 2.2, 3.3);

    initPhi(4.4, 5.5);
    initThetaPhi(6.6, 7.7);
    initEnergy(8.8, 9.9);

    Double_t x[3]={1.2, 2.3, 4.3};
    Double_t y[3]={1.2, 2.3, 4.3};
    Double_t z[3]={1.2, 2.3, 4.3};
    initXyz(x,y,z,3);

}

Bool_t AdfRecon::CompareInRange( const AdfRecon &ref, const std::string& name ) const {
    bool result = true;

    result = rootdatautil::CompareInRange(getEventNumber(),ref.getEventNumber(),"EventNumber") && result;
    result = rootdatautil::CompareInRange(getSpillNumber(),ref.getSpillNumber(),"SpillNumber") && result;
    result = rootdatautil::CompareInRange(getNumHighestCluster(),ref.getNumHighestCluster(),"numHighestsCluster") && result;
    result = rootdatautil::CompareInRange(getPx(),ref.getPx(),"Px") && result;
    result = rootdatautil::CompareInRange(getPy(),ref.getPy(),"Py") && result;
    result = rootdatautil::CompareInRange(getPz(),ref.getPz(),"Pz") && result;
    result = rootdatautil::CompareInRange(getPhiIn(),ref.getPhiIn(),"PhiIn") && result;
    result = rootdatautil::CompareInRange(getPhiOut(),ref.getPhiOut(),"PhiOut") && result;
    result = rootdatautil::CompareInRange(getTheta(),ref.getTheta(),"Theta") && result;
    result = rootdatautil::CompareInRange(getDeltaPhi(),ref.getDeltaPhi(),"DeltaPhi") && result;
    result = rootdatautil::CompareInRange(getReconstructedEnergy(),ref.getReconstructedEnergy(),"ReconstructedEnergy") && result;
    result = rootdatautil::CompareInRange(getCorrectedEnergy(),ref.getCorrectedEnergy(),"CorrectedEnergy") && result;

    result = rootdatautil::TObjArrayCompareInRange<reconRootData::TaggerCluster>(m_taggerClusterCol,ref.m_taggerClusterCol) && result;
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

} // end namespace
