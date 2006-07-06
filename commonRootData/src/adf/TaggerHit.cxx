
#include "commonRootData/adf/TaggerHit.h"
#include "commonRootData/RootDataUtil.h"
#include "Riostream.h"
#include "TCollection.h"  // Declares TIter

using namespace std;

ClassImp(commonRootData::TaggerHit)

namespace commonRootData {

TaggerHit::TaggerHit() {
    Clear();
}

TaggerHit::TaggerHit(UInt_t moduleId, UInt_t layerId, UInt_t stripId,
                     UInt_t pulseHght, Bool_t isPedSubtracted) {

    initialize(moduleId, layerId, stripId, pulseHght, isPedSubtracted);
}

TaggerHit::TaggerHit(const TaggerHit& copy):TObject(copy) {

    initialize(copy.m_moduleId, copy.m_layerId, copy.m_stripId, 
               copy.m_pulseHeight, copy.m_isPedestalSubtracted);
}


TaggerHit::~TaggerHit() {
    Clear();  
}

void TaggerHit::initialize(UInt_t moduleId, UInt_t layerId, UInt_t stripId,
                           UInt_t pulseHgt, Bool_t isPedSubtracted) {

    m_moduleId = moduleId;
    m_layerId = layerId;
    m_stripId = stripId;
    m_pulseHeight = pulseHgt;
    m_isPedestalSubtracted = isPedSubtracted;
}

TaggerHit& TaggerHit::operator=(const TaggerHit& copy) {
    initialize(copy.m_moduleId, copy.m_layerId, copy.m_stripId, 
               copy.m_pulseHeight, m_isPedestalSubtracted);
    return *this;
}


void TaggerHit::Clear(Option_t *option) {
    m_moduleId = 0;
    m_layerId = 0;
    m_stripId = 0;
    m_pulseHeight = 0;
    m_isPedestalSubtracted = false;
}

void TaggerHit::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "TaggerHit:" << endl;
    cout << "ModuleId: " << m_moduleId << " LayerId: " << m_layerId
         << " StripId: " << m_stripId << " PulseHeight: " << m_pulseHeight
         << " isPedestalSubtracted: " << m_isPedestalSubtracted << endl;
    cout << dec;
}

Bool_t TaggerHit::CompareInRange(const TaggerHit &ref, const std::string& name ) const{
    bool result = true;

    result = rootdatautil::CompareInRange(getModuleId(),ref.getModuleId(),"ModuleId") && result;
    result = rootdatautil::CompareInRange(getLayerId(),ref.getLayerId(),"LayerId") && result;
    result = rootdatautil::CompareInRange(getStripId(),ref.getStripId(),"StripId") && result;
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
