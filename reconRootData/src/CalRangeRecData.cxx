
// Include files
// Root Includes
#include "reconRootData/CalRangeRecData.h"

ClassImp(CalRangeRecData)
/*!
//------------------------------------------------------------------------------
//
// \class   CalRangeRecdata        
//  
// \brief   position, reconstructed energies for both faces of Xtal for Cal                                
//              
// Author:  A.Chekhtman, Apr, 17, 2002
//
//------------------------------------------------------------------------------
*/

// retrieve energy from specified face
Double_t CalRangeRecData::getEnergy(CalXtalId::XtalFace face) const {return face == CalXtalId::POS ? m_eneP : m_eneM;};

// retrieve energy range from specified face
Char_t CalRangeRecData::getRange(CalXtalId::XtalFace face) const {return face == CalXtalId::POS ? m_rangeP : m_rangeM;};
