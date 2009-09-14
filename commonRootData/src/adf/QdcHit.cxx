#include "commonRootData/adf/QdcHit.h"
#include "commonRootData/RootDataUtil.h"
#include "Riostream.h"
#include "TCollection.h"  // Declares TIter

ClassImp(commonRootData::QdcHit)

using namespace std;

namespace commonRootData {

QdcHit::QdcHit() {
    Clear();
}

QdcHit::QdcHit(UInt_t channel, Double_t pulseHght, UInt_t module,
               Double_t sigma, Bool_t isPedSubtracted) {

    initialize(channel, pulseHght, module, sigma, isPedSubtracted);
}

QdcHit::QdcHit(const QdcHit& copy):TObject(copy) {

    initialize(copy.m_channel, copy.m_pulseHeight, copy.m_module, 
               copy.m_sigma, copy.m_isPedestalSubtracted);
}


QdcHit::~QdcHit() {
    Clear();  
}

void QdcHit::initialize(UInt_t channel, Double_t pulseHgt, UInt_t module,
                        Double_t sigma, Bool_t isPedSubtracted) {

    m_channel = channel;
    m_pulseHeight = pulseHgt;
    m_sigma = sigma;
    m_module = module;
    m_isPedestalSubtracted = isPedSubtracted;
}

QdcHit& QdcHit::operator=(const QdcHit& copy) {
    initialize(copy.m_channel, copy.m_pulseHeight, copy.m_module, 
               copy.m_sigma,copy.m_isPedestalSubtracted);
    return *this;
}


void QdcHit::Clear(Option_t *option) {
    TObject::Clear(option);
    m_channel = 0;
    m_module = 0;
    m_pulseHeight = 0.0;
    m_sigma = 0.0;
    m_isPedestalSubtracted = false;
}

void QdcHit::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "QdcHit:" << endl;
    cout << "Channel: " << m_channel  
         << " PulseHeight: " << m_pulseHeight
         << " Sigma: " << m_sigma
         << " Module: " << m_module
         << " isPedestalSubtracted: " << m_isPedestalSubtracted << endl;
    cout << dec;
}

Bool_t QdcHit::CompareInRange( const QdcHit &ref, const std::string& name ) const {
    bool result = true;

    result = rootdatautil::CompareInRange(getChannel(),ref.getChannel(),"Channel") && result;
    result = rootdatautil::CompareInRange(getPulseHeight(),ref.getPulseHeight(),"PulseHeight") && result;
    result = rootdatautil::CompareInRange(isPedestalSubtracted(),ref.isPedestalSubtracted(),"IsPedestalSubtracted") && result;
    result = rootdatautil::CompareInRange(getModule(), ref.getModule(),"Module") && result;
    result = rootdatautil::CompareInRange(getSigma(), ref.getSigma(),"Sigma") && result;

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
