#include "digiRootData/CalLogReadout.h"

ClassImp(CalLogReadout)

/*!
//------------------------------------------------------------------------------
//
// \class   CalLogReadout        
//  
// \brief   Pulse heights and energy range for both faces of log for Cal                                
//              
// Author:  J. Eric Grove, 23 Feb 2001
//
//------------------------------------------------------------------------------
*/
/// Destructor
CalLogReadout::~CalLogReadout() {}

// retrieve pulse height from specified face
Short_t CalLogReadout::getAdc(CalLogReadout::LogFace face) const {return face == POS ? m_adcP : m_adcM;}

// retrieve energy range from specified face
Char_t CalLogReadout::getRange(CalLogReadout::LogFace face) const {return face == POS ? m_rangeP : m_rangeM;}
