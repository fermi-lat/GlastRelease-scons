#include "digiRootData/CalXtalReadout.h"
#include <iostream>

ClassImp(CalXtalReadout)

// Class   CalXtalReadout        
// Description Pulse heights and energy range for both faces of log for Cal.
// Based on the CalXtalReadout TDS class defined in GlastEvent 
// written by J. Eric Grove.
//              
// Author: Heather Kelly
//

CalXtalReadout::CalXtalReadout() {
    Clear();
}

CalXtalReadout::CalXtalReadout(Char_t rangeP, UShort_t adcP, Char_t rangeM, 
                               UShort_t adcM) :
m_adcP(adcP), m_adcM(adcM), m_rangeP(rangeP), m_rangeM(rangeM)
{
}

CalXtalReadout::~CalXtalReadout() {
}

void CalXtalReadout::initialize(Char_t rangeP, UShort_t adcP, Char_t rangeM, 
                          UShort_t adcM) {
    Clear();
    m_rangeP = rangeP;
    m_adcP = adcP;
    m_rangeM = rangeM;
    m_adcM = adcM;
}

void CalXtalReadout::Clear(Option_t *option) {
    m_rangeP = 0;
    m_adcP = 0;
    m_rangeM = 0;
    m_adcM = 0;
}

void CalXtalReadout::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "POS range, adc " << Short_t(m_rangeP) << " " << m_adcP << endl;
    cout << "NEG range, adc " << Short_t(m_rangeM) << " " << m_adcM << endl;
}


UShort_t CalXtalReadout::getAdc(CalXtalId::XtalFace face) const {
    return face == CalXtalId::POS ? m_adcP : m_adcM;
}

Char_t CalXtalReadout::getRange(CalXtalId::XtalFace face) const {
    return face == CalXtalId::POS ? m_rangeP : m_rangeM;
}
