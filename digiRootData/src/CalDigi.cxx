#include "digiRootData/CalDigi.h"
ClassImp(CalDigi)

//
// Class:  CalDigi        
// Description: Digitizations for Cal
// Based on the CalDigi TDS class in GlastEvent written by J. Eric Grove.
//              
// Author:  Heather Kelly
//

CalDigi::CalDigi() {
    m_readoutCol = 0;
    Clear();
}

CalDigi::~CalDigi() { 
    Clear();
    if(m_readoutCol) {
        m_readoutCol->Delete();
        delete m_readoutCol;
        m_readoutCol = 0;
    }
}

void CalDigi::initialize(CalXtalId::CalTrigMode m, const CalXtalId &id) {
    m_mode = m;
    m_xtalId = id;
    m_readoutCol = new TClonesArray("CalXtalReadout", 1);
}

void CalDigi::Clear(Option_t *option) {
    m_numXtals = 0;
    m_xtalId.Clear();
    if (m_readoutCol) m_readoutCol->Clear();
    m_mode = CalXtalId::BESTRANGE;
}

void CalDigi::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Id: " << endl;
    m_xtalId.Print();
    cout << "Num Readouts: " << m_readoutCol->GetEntries() << endl;
}

const CalXtalReadout* CalDigi::addReadout(Char_t rangeP, UShort_t adcP, Char_t rangeM, UShort_t adcM) 
{ 
    // Add a new CalXtalReadout entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    //++m_numXtals;
    TClonesArray &xtals = *m_readoutCol;
    new(xtals[m_numXtals]) CalXtalReadout(rangeP, adcP, rangeM, adcM);
    return ((CalXtalReadout*)(xtals[m_numXtals++]));
} 

Char_t CalDigi::getRange(short readoutIndex, CalXtalId::XtalFace face) const
{
    return (readoutIndex < m_readoutCol->GetEntries()) 
        ? ((CalXtalReadout*)m_readoutCol->At(readoutIndex))->getRange(face) : (Char_t)-1;
}

Short_t CalDigi::getAdc(Short_t readoutIndex, CalXtalId::XtalFace face) const
{
    return (readoutIndex < m_readoutCol->GetEntries()) 
        ? ((CalXtalReadout*)m_readoutCol->At(readoutIndex))->getAdc(face) : (Short_t)-1;
}

const CalXtalReadout* CalDigi::getXtalReadout(Short_t readoutIndex)
{
    if ( readoutIndex < m_readoutCol->GetEntries() )
        return ((CalXtalReadout*)m_readoutCol->At(readoutIndex));
    else
        return 0;
    
}

Short_t CalDigi::getAdcSelectedRange(Char_t range, CalXtalId::XtalFace face) const
{
    Char_t nRanges = (Char_t)m_readoutCol->GetEntries();
    if (nRanges <= 0) return (Short_t)-1;
    CalXtalReadout* firstXtal = (CalXtalReadout*)m_readoutCol->At(0);

    if (nRanges == 1) {
        return (range == firstXtal->getRange(face)) 
        ? firstXtal->getAdc(face) : (Short_t)-1;
    } else {
        UInt_t index = (nRanges + range - firstXtal->getRange(face)) % nRanges;
        CalXtalReadout* xtal = (CalXtalReadout*) m_readoutCol->At(index);
        return xtal->getAdc(face);
    }
}
