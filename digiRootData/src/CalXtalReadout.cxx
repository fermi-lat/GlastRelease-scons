#include "digiRootData/CalXtalReadout.h"

ClassImp(CalXtalReadout)

// Class   CalXtalReadout        
// Description Pulse heights and energy range for both faces of log for Cal   
//              
// Author:  J. Eric Grove, 23 Feb 2001
//

CalXtalReadout::CalXtalReadout() {
    Clear();
}

CalXtalReadout::CalXtalReadout(Char_t rangeP, UShort_t adcP, Char_t rangeM, 
                               UShort_t adcM) :
m_rangeP(rangeP), m_adcP(adcP), m_rangeM(adcM), m_adcM(adcM)
{
}

CalXtalReadout::~CalXtalReadout() {
}

void CalXtalReadout::init(Char_t rangeP, UShort_t adcP, Char_t rangeM, 
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

}


Short_t CalXtalReadout::getAdc(CalXtalId::XtalFace face) const {
    return face == CalXtalId::POS ? m_adcP : m_adcM;
}

Char_t CalXtalReadout::getRange(CalXtalId::XtalFace face) const {
    return face == CalXtalId::POS ? m_rangeP : m_rangeM;
}
