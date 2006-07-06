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

QdcHit::QdcHit(UInt_t channel, UInt_t pulseHght, Bool_t isPedSubtracted) {

    initialize(channel, pulseHght, isPedSubtracted);
}

QdcHit::QdcHit(const QdcHit& copy):TObject(copy) {

    initialize(copy.m_channel, copy.m_pulseHeight, copy.m_isPedestalSubtracted);
}


QdcHit::~QdcHit() {
    Clear();  
}

void QdcHit::initialize(UInt_t channel, UInt_t pulseHgt, Bool_t isPedSubtracted) {

    m_channel = channel;
    m_pulseHeight = pulseHgt;
    m_isPedestalSubtracted = isPedSubtracted;
}

QdcHit& QdcHit::operator=(const QdcHit& copy) {
    initialize(copy.m_channel, m_pulseHeight, copy.m_isPedestalSubtracted);
    return *this;
}


void QdcHit::Clear(Option_t *option) {
    m_channel = 0;
    m_pulseHeight = 0;
    m_isPedestalSubtracted = false;
}

void QdcHit::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "QdcHit:" << endl;
    cout << "Channel: " << m_channel  
         << " PulseHeight: " << m_pulseHeight
         << " isPedestalSubtracted: " << m_isPedestalSubtracted << endl;
    cout << dec;
}

Bool_t QdcHit::CompareInRange( const QdcHit &ref, const std::string& name ) const {
    bool result = true;

    result = rootdatautil::CompareInRange(getChannel(),ref.getChannel(),"Channel") && result;
    result = rootdatautil::CompareInRange(getPulseHeight(),ref.getPulseHeight(),"PulseHeight") && result;
    result = rootdatautil::CompareInRange(isPedestalSubtracted(),ref.isPedestalSubtracted(),"IsPedestalSubtracted") && result;

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
