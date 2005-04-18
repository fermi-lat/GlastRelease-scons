
#include "reconRootData/CalRangeRecData.h"
#include <iostream>

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
