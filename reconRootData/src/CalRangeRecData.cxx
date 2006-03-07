
#include <reconRootData/CalRangeRecData.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(CalRangeRecData)

//
// CalRangeRecData        
//  
// position and reconstructed energies for both faces of Xtal for Cal                                
//              
// Author:  A.Chekhtman, Apr, 17, 2002
//

// retrieve energy from specified face
Double_t CalRangeRecData::getEnergy(CalXtalId::XtalFace face) const {return face == CalXtalId::POS ? m_eneP : m_eneM;};

// retrieve energy range from specified face
Char_t CalRangeRecData::getRange(CalXtalId::XtalFace face) const {return face == CalXtalId::POS ? m_rangeP : m_rangeM;};

void CalRangeRecData::Clear(Option_t* /* option */) {
    m_eneP = 0.0;
    m_eneM = 0.0;
    m_pos = TVector3(0.0, 0.0, 0.0);
    m_rangeP = 0;
    m_rangeM = 0;
}

void CalRangeRecData::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout << "POS (range, energy): (" << Int_t(m_rangeP) << "," << m_eneP << ")" << endl;
    cout << "NEG (range, energy): (" << Int_t(m_rangeM) << "," << m_eneM << ")" << endl;
}

// dummy data, just for tests
void CalRangeRecData::Fake( Int_t ievent, UInt_t /*rank*/, Float_t randNum ) {
    
    Float_t f = Float_t(ievent) ;
    m_pos = TVector3(4.5,7.5,8.5) ;
    m_rangeP = CalXtalId::LEX8 ;
    m_eneP = randNum*f ;
    m_rangeM = CalXtalId::HEX8 ;
    m_eneM = randNum*4.0 ;

 }

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att,ref.get ## att,#att)

Bool_t CalRangeRecData::CompareInRange( const CalRangeRecData & ref, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_IN_RANGE(Position()) && result ;
    result = COMPARE_IN_RANGE(Energy(CalXtalId::POS)) && result ;
    result = COMPARE_IN_RANGE(Energy(CalXtalId::NEG)) && result ;
    result = COMPARE_IN_RANGE(Range(CalXtalId::POS)) && result ;
    result = COMPARE_IN_RANGE(Range(CalXtalId::NEG)) && result ;

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
