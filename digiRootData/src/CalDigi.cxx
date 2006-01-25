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
//    m_readoutCol = 0;
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
    for(int i=0;i<4;i++)(m_readoutArr[i]).Clear();
    //if (m_readoutCol) m_readoutCol->Clear();
    m_mode = CalXtalId::BESTRANGE;
}

void CalDigi::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Id: " << endl;
    m_xtalId.Print();
    // Maintain Backward Compatibility
    //if (m_readoutCol)
    //    cout << "Num Readouts: " << m_readoutCol->GetEntries() << endl;
    //else
        cout << "Num Readouts: " << m_numReadouts << endl;
}

const CalXtalReadout* CalDigi::addReadout(Char_t rangeP, UShort_t adcP, Char_t rangeM, UShort_t adcM) 
{ 
    // Purpose and Method: Add a new CalXtalReadout entry
    //++m_numXtals;
    m_readoutArr[m_numReadouts] = CalXtalReadout(rangeP, adcP, rangeM, adcM);
    return &(m_readoutArr[m_numReadouts++]);
} 

Char_t CalDigi::getRange(UShort_t readoutIndex, CalXtalId::XtalFace face) const
{
    // Maintain backward compatibility
    //if (m_readoutCol) {
    //    return (readoutIndex < m_readoutCol->GetEntries())
    //        ? ((CalXtalReadout*)m_readoutCol->At(readoutIndex))->getRange(face) : (Char_t)-1;
   // }

    return (readoutIndex < m_numReadouts) 
        ? (m_readoutArr[readoutIndex]).getRange(face) : (Char_t)-1;
}

Short_t CalDigi::getAdc(UShort_t readoutIndex, CalXtalId::XtalFace face) const
{
    // Maintain backward compatibility
    //if (m_readoutCol) {
    //    return (readoutIndex < m_readoutCol->GetEntries())
    //        ? ((CalXtalReadout*)m_readoutCol->At(readoutIndex))->getAdc(face) : (Short_t)-1;
    //}

    return (readoutIndex < m_numReadouts) 
        ? (m_readoutArr[readoutIndex]).getAdc(face) : (Short_t)-1;
}

const CalXtalReadout* CalDigi::getXtalReadout(UShort_t readoutIndex) const
{
    // Maintain Backward Compatibility
    //if (m_readoutCol) {
    //    return (CalXtalReadout*)m_readoutCol->At(readoutIndex);
    //}
    if ( readoutIndex < m_numReadouts )
        return &(m_readoutArr[readoutIndex]);
    else
        return 0;
    
}

Short_t CalDigi::getAdcSelectedRange(Char_t range, CalXtalId::XtalFace face) const
{
    // Maintain Backward Compatibility
    //if (m_readoutCol) {
    //    if (m_readoutCol->GetEntries() <= 0) return (Short_t)-1;
    //    CalXtalReadout* firstXtal = (CalXtalReadout*)m_readoutCol->At(0);
        
    //    Char_t nRanges = (Char_t)m_readoutCol->GetEntries();
    //    if (nRanges == 1) {
    //        return (range == firstXtal->getRange(face)) 
    //            ? firstXtal->getAdc(face) : (Short_t)-1;
    //    } else {
            //UInt_t index = (nRanges + range - firstXtal->getRange(face)) % nRanges;
            //CalXtalReadout* xtal = (CalXtalReadout*) m_readoutCol->At(index);
            //return xtal->getAdc(face);
      //  }

    //}

    if (m_numReadouts <= 0) return (Short_t)-1;
    const CalXtalReadout* firstXtal = &(m_readoutArr[0]);

    if (m_numReadouts == 1) {
        return (range == firstXtal->getRange(face)) 
        ? firstXtal->getAdc(face) : (Short_t)-1;
    } else {
        UInt_t index = (m_numReadouts + range - firstXtal->getRange(face)) % m_numReadouts;
        const CalXtalReadout* xtal = & (m_readoutArr[index]);
        return xtal->getAdc(face);
    }
}


const CalXtalReadout* CalDigi::getReadoutCol() const { 
    // Purpose and Method:  Provide access to the full array of CalXtalReadouts
    // For backward compatibility - we handle the possibility that the TClonesArray* version
    // is available - in which case we fill the fixed array using the quantities in the 
    // m_readoutCol (TClonesArray) data member.  This conversion should happen just the first time
    // that getReadoutCol is called for a particular event.
    using namespace std;
    //if ((m_readoutCol) && (!m_numReadouts)){
    //    TIter cloneIter(m_readoutCol);
    //    CalXtalReadout *readout = 0;
    //    while ((readout = (CalXtalReadout*)cloneIter.Next())!=0) {
    //        addReadout(readout->getRange(CalXtalId::POS), readout->getAdc(CalXtalId::POS),
    //            readout->getRange(CalXtalId::NEG), readout->getAdc(CalXtalId::NEG));
    //    }

    //}
        
    return m_readoutArr; 
}
