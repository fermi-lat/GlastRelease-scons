
#include "reconRootData/CalXtalRecData.h"
#include <commonRootData/RootDataUtil.h>
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

void CalXtalRecData::Clear(Option_t* /* option */) {
    m_recData.clear();
}

void CalXtalRecData::Print(Option_t* /* option */) const {
    using namespace std;
    TObject::Print();
    cout << "XtalId (tower, layer, column): " << m_xtalId.getTower() << "," 
        << m_xtalId.getLayer() << ", " << m_xtalId.getColumn() << endl;
    cout << "BESTRANGE POS(range, energy): (" << Int_t(getRange(0, CalXtalId::POS))
        << "," << getEnergy(0, CalXtalId::POS) << ")"
        << "  NEG (range, energy): (" << Int_t(getRange(0, CalXtalId::NEG))
        << "," << getEnergy(0, CalXtalId::NEG) << ")" << std::endl;

}


Char_t CalXtalRecData::getRange(UShort_t readoutIndex, CalXtalId::XtalFace face) const 
{
    return (readoutIndex < m_recData.size()) ? ((m_recData[readoutIndex])).getRange(face) : (char)-1;
}


Double_t CalXtalRecData::getEnergy(UShort_t readoutIndex, CalXtalId::XtalFace face) const
{
    return (readoutIndex < m_recData.size()) ? ((m_recData[readoutIndex])).getEnergy(face) : -1.0;
};


Double_t CalXtalRecData::getEnergy() const
{
    return (getEnergy(0,CalXtalId::POS)
        +getEnergy(0,CalXtalId::NEG))/2;
};


const CalRangeRecData* CalXtalRecData::getRangeRecData(UShort_t readoutIndex) const
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


//===================================================
// dummy data, just for tests
//===================================================

void CalXtalRecData::Fake( Int_t ievent, UInt_t rank, Float_t randNum ) {

    CalXtalId id;
    id.init(1, 2, 3);
    initialize(CalXtalId::BESTRANGE,id) ;
    CalRangeRecData rec ;
    rec.Fake(ievent,rank,randNum) ;
    addRangeRecData(rec) ;

 }

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att,ref.get ## att,#att)

Bool_t CalXtalRecData::CompareInRange( const CalXtalRecData & ref, const std::string & name ) const {

    bool result = true ;

    // comparison with ref
    
    result = COMPARE_IN_RANGE(Mode()) && result ;
    result = COMPARE_IN_RANGE(PackedId()) && result ;
    result = COMPARE_IN_RANGE(Position()) && result ;

    result = COMPARE_IN_RANGE(Range(0,CalXtalId::POS)) && result ;
    result = COMPARE_IN_RANGE(Range(0,CalXtalId::NEG)) && result ;

    result = COMPARE_IN_RANGE(Energy(0,CalXtalId::POS)) && result ;
    result = COMPARE_IN_RANGE(Energy(0,CalXtalId::NEG)) && result ;
    
    result = getRangeRecData(0)->CompareInRange(*ref.getRangeRecData(0)) && result ;

    // additionnal consistency checks
    
    Double_t energy = getEnergy();
    Double_t energyP = getEnergy(0,CalXtalId::POS) ;
    Double_t energyN = getEnergy(0,CalXtalId::NEG) ;
    result = result && rootdatautil::CompareInRange(energy,
      (energyP+energyN)/2.,"Average Energy") ;
      
    const CalRangeRecData * rangeRecData = getRangeRecData(0) ;
    result = result && rootdatautil::CompareInRange(energyP,
      rangeRecData->getEnergy(CalXtalId::POS),"CalRangeRecData POS Energy") ;
    result = result && rootdatautil::CompareInRange(energyN,
      rangeRecData->getEnergy(CalXtalId::NEG),"CalRangeRecData NEG Energy") ;
      
    result = result && rootdatautil::CompareInRange(energyP,
      getEnergySelectedRange(CalXtalId::LEX8,CalXtalId::POS),"Selected POS Energy") ;
    result = result && rootdatautil::CompareInRange(energyN,
      getEnergySelectedRange(CalXtalId::HEX8,CalXtalId::NEG),"Selected NEG Energy") ;

    // check range and faces that should return -1, since they are undefined
    result = result && rootdatautil::CompareInRange(-1.,
      getEnergySelectedRange(CalXtalId::LEX1,CalXtalId::POS),"Undefined Range Face Pair") ;
    result = result && rootdatautil::CompareInRange(-1.,
      getEnergySelectedRange(CalXtalId::HEX1,CalXtalId::POS),"Undefined Range Face Pair") ;
    result = result && rootdatautil::CompareInRange(-1.,
      getEnergySelectedRange(CalXtalId::HEX8,CalXtalId::POS),"Undefined Range Face Pair") ;
    result = result && rootdatautil::CompareInRange(-1.,
      getEnergySelectedRange(CalXtalId::LEX1,CalXtalId::NEG),"Undefined Range Face Pair") ;
    result = result && rootdatautil::CompareInRange(-1.,
      getEnergySelectedRange(CalXtalId::LEX8,CalXtalId::NEG),"Undefined Range Face Pair") ;
    result = result && rootdatautil::CompareInRange(-1.,
      getEnergySelectedRange(CalXtalId::HEX1,CalXtalId::NEG),"Undefined Range Face Pair") ;

    // end
    
    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}

