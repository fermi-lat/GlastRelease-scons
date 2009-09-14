
#include "commonRootData/adf/ScalerHit.h"
#include "commonRootData/RootDataUtil.h"
#include "Riostream.h"
#include "TCollection.h"  // Declares TIter

using namespace std;

ClassImp(commonRootData::ScalerHit)

namespace commonRootData {

ScalerHit::ScalerHit() {
    Clear();
}

ScalerHit::ScalerHit(UInt_t channel, UInt_t value) {

    initialize(channel, value);
}

ScalerHit::ScalerHit(const ScalerHit& copy):TObject(copy) {

    initialize(copy.m_channel, copy.m_value);
}


ScalerHit::~ScalerHit() {
    Clear();  
}

void ScalerHit::initialize(UInt_t channel, UInt_t value) {

    m_channel = channel;
    m_value = value;
}

ScalerHit& ScalerHit::operator=(const ScalerHit& copy) {
    initialize(copy.m_channel, copy.m_value);  
    return *this;
}


void ScalerHit::Clear(Option_t *option) {
    TObject::Clear(option);
    m_channel = 0;
    m_value = 0;
}

void ScalerHit::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "ScalerHit:" << endl;
    cout << "Channel: " << m_channel << " Value: " << m_value << endl; 
    cout << dec;
}

Bool_t ScalerHit::CompareInRange(const ScalerHit &ref, const std::string& name ) const{
    bool result = true;

    result = rootdatautil::CompareInRange(getChannel(),ref.getChannel(),"Channel") && result;
    result = rootdatautil::CompareInRange(getValue(),ref.getValue(),"Value") && result;

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
