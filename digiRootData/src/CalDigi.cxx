#include "digiRootData/CalDigi.h"
#include <iostream>

ClassImp(CalDigi)

//
// Class:  CalDigi        
// Description: Digitizations for Cal
// Based on the CalDigi TDS class in GlastEvent written by J. Eric Grove.
//              
// Author:  Heather Kelly
//

CalDigi::CalDigi() {
    Clear();
}

CalDigi::~CalDigi() { 
    Clear();    
}

void CalDigi::initialize(CalXtalId::CalTrigMode m, const CalXtalId &id) {
    m_mode = m;
    m_xtalId = id;
}

void CalDigi::Clear(Option_t *option) {
    m_numReadouts = 0;
    m_xtalId.Clear();
    for(int i=0;i<4;i++)(m_readoutCol[i]).Clear();
    m_mode = CalXtalId::BESTRANGE;
}

void CalDigi::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Id: " << endl;
    m_xtalId.Print();
    cout << "Num Readouts: " << m_numReadouts << endl;
}

const CalXtalReadout* CalDigi::addReadout(Char_t rangeP, UShort_t adcP, Char_t rangeM, UShort_t adcM) 
{ 
    // Add a new CalXtalReadout entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    //++m_numXtals;
    m_readoutCol[m_numReadouts] =  CalXtalReadout(rangeP, adcP, rangeM, adcM);
    return &(m_readoutCol[m_numReadouts++]);
} 

Char_t CalDigi::getRange(short readoutIndex, CalXtalId::XtalFace face) const
{
    return (readoutIndex < m_numReadouts) 
        ? (m_readoutCol[readoutIndex]).getRange(face) : (Char_t)-1;
}

Short_t CalDigi::getAdc(Short_t readoutIndex, CalXtalId::XtalFace face) const
{
    return (readoutIndex < m_numReadouts) 
        ? (m_readoutCol[readoutIndex]).getAdc(face) : (Short_t)-1;
}

const CalXtalReadout* CalDigi::getXtalReadout(Short_t readoutIndex)
{
    if ( readoutIndex < m_numReadouts )
        return &(m_readoutCol[readoutIndex]);
    else
        return 0;
    
}

Short_t CalDigi::getAdcSelectedRange(Char_t range, CalXtalId::XtalFace face) const
{
    if (m_numReadouts <= 0) return (Short_t)-1;
    const CalXtalReadout* firstXtal = &(m_readoutCol[0]);

    if (m_numReadouts == 1) {
        return (range == firstXtal->getRange(face)) 
        ? firstXtal->getAdc(face) : (Short_t)-1;
    } else {
        UInt_t index = (m_numReadouts + range - firstXtal->getRange(face)) % m_numReadouts;
        const CalXtalReadout* xtal = & (m_readoutCol[index]);
        return xtal->getAdc(face);
    }
}
