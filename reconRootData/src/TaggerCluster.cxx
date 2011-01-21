
#include "reconRootData/TaggerCluster.h"
#include "commonRootData/RootDataUtil.h"
#include "Riostream.h"
#include "TCollection.h"  // Declares TIter

using namespace std;

ClassImp(reconRootData::TaggerCluster)

namespace reconRootData {

TaggerCluster::TaggerCluster() {
    m_taggerHitCol = 0; 
    m_numTaggerHit = -1;

    m_baricenterPosition = 0.0;
    m_totalPulseHeight = 0.0;
    m_eta = 0.0;
    m_highestHit.Clear();
    m_x = 0.0;
    m_y = 0.0;
    m_z = 0.0;
}


TaggerCluster::~TaggerCluster() {
  
        if (m_taggerHitCol) {
        m_taggerHitCol->Delete();
        delete m_taggerHitCol;
        m_taggerHitCol = 0;
        }
  
}


void TaggerCluster::initialize(Double_t pos, Double_t pulseHeight, Double_t eta
) {
    m_baricenterPosition = pos;
    m_totalPulseHeight = pulseHeight;
    m_eta = eta;
}

void TaggerCluster::initXyz(Double_t x, Double_t y, Double_t z) {
    m_x = x;
    m_y = y;
    m_z = z; 
}

void TaggerCluster::Clear(Option_t *option) {
    TObject::Clear(option);

    if (m_taggerHitCol) m_taggerHitCol->Clear("C");
    m_numTaggerHit = -1;

    m_baricenterPosition = 0.0;
    m_totalPulseHeight = 0.0;
    m_eta = 0.0;
    m_highestHit.Clear();
    m_x = 0.0;
    m_y = 0.0;
    m_z = 0.0;
}

void TaggerCluster::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "TaggerCluster:" << endl;
    cout << "Pos: " << m_baricenterPosition << " PulseHeight: " 
         << m_totalPulseHeight << " Eta: " << m_eta  
         << "(x,y,z): (" << m_x << "," << m_y << "," << m_z << ")" << endl;
    cout << "Num Tagger Hits: " << m_taggerHitCol->GetEntriesFast() <<endl;
    cout << dec;
}


const commonRootData::TaggerHit* TaggerCluster::getTaggerHit(UInt_t ind) const {
    if (((Int_t)ind) < m_taggerHitCol->GetEntriesFast()) 
        return ((commonRootData::TaggerHit*)m_taggerHitCol->At(ind));
    else
        return 0;
}

commonRootData::TaggerHit* TaggerCluster::addTaggerHit(UInt_t moduleId, 
      UInt_t layerId, UInt_t stripId, Double_t pulseHgt, Double_t sigma,
      Bool_t isPedSub) {
    // Add a new TaggerHit entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    if (!m_taggerHitCol) m_taggerHitCol = new TClonesArray("commonRootData::TaggerHit",1);
    ++m_numTaggerHit;
    TClonesArray &localCol = *m_taggerHitCol;
    new(localCol[m_numTaggerHit]) commonRootData::TaggerHit(moduleId, layerId, stripId, pulseHgt, sigma, isPedSub);
    return ((commonRootData::TaggerHit*)(localCol[m_numTaggerHit]));
}


void TaggerCluster::Fake( Int_t /*ievent*/, Float_t /*randNum*/ ) {
    initialize(1.1, 2.2, 3.3);
    initXyz(4.4, 5.5, 6.6);
    // Add 2 TaggerHits
    addTaggerHit(1,2,3,4.4,5.5,true);
    addTaggerHit(5,6,7,8.8,6.6,false);
   
}

Bool_t TaggerCluster::CompareInRange( const TaggerCluster &ref, const std::string& name ) const {
    bool result = true;

    result = rootdatautil::CompareInRange(getPosition(),ref.getPosition(),"Position") && result;
    result = rootdatautil::CompareInRange(getPulseHeight(),ref.getPulseHeight(),"PulseHeight") && result;
    result = rootdatautil::CompareInRange(getEta(),ref.getEta(),"Eta") && result;
    result = rootdatautil::CompareInRange(getX(),ref.getX(),"X") && result;
    result = rootdatautil::CompareInRange(getY(),ref.getY(),"Y") && result;
    result = rootdatautil::CompareInRange(getZ(),ref.getZ(),"Z") && result;

    result = rootdatautil::TObjArrayCompareInRange<commonRootData::TaggerHit>(m_taggerHitCol,ref.m_taggerHitCol) && result;

    if (!result) {
        if (name == "") {
            std::cout<<"Comparison ERROR for " << ClassName() << std::endl;
        } else {
            std::cout << "Comparison ERROR for " << name<<std::endl;
        }
    }

    return result ;

}

}
