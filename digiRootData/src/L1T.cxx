//                                                                       
// The L1T class contains the information about the Level 1 Trigger count
// for an event.  
//                                                                       

#include "digiRootData/L1T.h"
#include <iostream>

ClassImp(L1T)

L1T::L1T() : m_trigger(0)
{
}

L1T::L1T(UInt_t trigger) : m_trigger(trigger) {

}

L1T::L1T(const L1T& level1) {
    m_trigger = level1.m_trigger;
}

L1T::~L1T() {

}

void L1T::initialize(UInt_t trigger) {
    m_trigger = trigger;
}

void L1T::Clear(Option_t *option) {
    m_trigger = 0;
}

void L1T::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Trigger Bits: " << m_trigger << endl;
}
