
#include "reconRootData/CalXtalRecData.h"
#include <iostream>

ClassImp(CalXtalRecData)

// CalXtalRecData        
// reconstructed data for a calorimeter crystal                                
//
// Author: A.Chekhtman


void CalXtalRecData::initialize(CalXtalId::CalTrigMode m, CalXtalId id) {
    m_mode = m;
    m_xtalId = id;
}

void CalXtalRecData::Clear(Option_t *option) {
    m_recData.clear();
}

void CalXtalRecData::Print(Option_t *option) {
    using namespace std;
    TObject::Print();
    cout << "XtalId (tower, layer, column): " << m_xtalId.getTower() << "," 
        << m_xtalId.getLayer() << ", " << m_xtalId.getColumn() << endl;

}


Char_t CalXtalRecData::getRange(Short_t readoutIndex, CalXtalId::XtalFace face) const 
{
    return (readoutIndex < m_recData.size()) ? ((m_recData[readoutIndex])).getRange(face) : (char)-1;
}


Double_t CalXtalRecData::getEnergy(Short_t readoutIndex, CalXtalId::XtalFace face) const
{
    return (readoutIndex < m_recData.size()) ? ((m_recData[readoutIndex])).getEnergy(face) : -1.0;
};


Double_t CalXtalRecData::getEnergy()
{
    return (getEnergy(0,CalXtalId::POS)
        +getEnergy(0,CalXtalId::NEG))/2;
};


CalRangeRecData* CalXtalRecData::getRangeRecData(Short_t readoutIndex)
{
    if ( readoutIndex < m_recData.size() )
        return &(m_recData[readoutIndex]);
    else
        return 0;
    
};

Double_t CalXtalRecData::getEnergySelectedRange(CalXtalId::AdcRange range, CalXtalId::XtalFace face) const
{
    Int_t nRanges = m_recData.size();
    if (nRanges == 1)
        return (range == ((m_recData[0])).getRange(face)) ? ((m_recData[0])).getEnergy(face) : -1.0;
    else
        return ((m_recData[(nRanges + range - ((m_recData[0])).getRange(face)) % nRanges])).getEnergy(face);
};

