//                                                                       
// The L1T class contains the information about the Level 1 Trigger count
// for an event.  
//                                                                       

#include "digiRootData/L1T.h"
#include <iostream>

ClassImp(L1T)

L1T::L1T() : m_trigger(0)
{
    Clear();
}

L1T::L1T(UInt_t trigger, UInt_t *triRowBits) {
    Clear();
    initialize(trigger, triRowBits);
}

L1T::L1T(const L1T& level1):TObject(level1) {
    Clear();
    m_trigger = level1.m_trigger;
    Int_t iTower;
    for (iTower = 0; iTower < 16; iTower++) 
        m_triRowBits[iTower] = level1.m_triRowBits[iTower];
}

L1T::~L1T() {

}

void L1T::initialize(UInt_t trigger, UInt_t *triRowBits) {
    m_trigger = trigger;
    Int_t tower;
    for (tower = 0; tower < 16; tower++) 
        setTriRowBits(tower, triRowBits[tower]); 
}

void L1T::setTriRowBits(UInt_t tower, UInt_t bits) {
    if ( (tower < 0) || (tower > 15)) return;
    m_triRowBits[tower] = bits;
}

void L1T::setTriRowBits(UInt_t *bits) {
    Int_t itower;
    for (itower = 0; itower < 16; itower++) 
        setTriRowBits(itower, bits[itower]);
}


void L1T::Clear(Option_t *option) {
    m_trigger = 0;
    Int_t iTower;
    for (iTower = 0; iTower < 16; iTower++) 
        m_triRowBits[iTower] = 0;
  
}

void L1T::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Trigger Bits: " << m_trigger << endl;
    cout << "TriRowBits {Tower 0 through 16} : " << endl;
    Int_t iTower;
    for (iTower = 0; iTower < 8; iTower++)
        cout << m_triRowBits[iTower] << " ";
    cout << endl;
    for (iTower = 8; iTower < 16; iTower++)
        cout << m_triRowBits[iTower] << " ";
    cout << endl;
   
}

UInt_t L1T::getTriRowBits(const Int_t tower) const {
    if ((tower < 0) || (tower > 15)) return 0;

    return m_triRowBits[tower];
}

