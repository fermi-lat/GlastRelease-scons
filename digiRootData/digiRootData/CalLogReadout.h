
#ifndef digiRootData_CalLogReadout_H
#define digiRootData_CalLogReadout_H 1


// Include files
#include <iostream>
#include <vector>
#include "TObject.h"


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
class CalLogReadout : public TObject  { 
    
public:
    
    CalLogReadout() {};
    CalLogReadout(Char_t rangeP, Short_t adcP, Char_t rangeM, Short_t adcM) :
    m_rangeP(rangeP), 
        m_adcP(adcP), 
        m_rangeM(adcM), 
        m_adcM(adcM)
    {};
    
    /// Destructor
    ~CalLogReadout();
    
    /// log ends are labeled by POSitive or NEGative face
    typedef enum
    {
        POS = 0,
        NEG
    } LogFace;
    
    // retrieve pulse height from specified face
    Short_t getAdc(CalLogReadout::LogFace face) const;
    
    // retrieve energy range from specified face
    Char_t getRange(CalLogReadout::LogFace face) const;
    
private:
    
    Short_t m_adcP;
    Short_t m_adcM;
    Char_t  m_rangeP;
    Char_t  m_rangeM;
    
    ClassDef(CalLogReadout,1)
};



#endif
