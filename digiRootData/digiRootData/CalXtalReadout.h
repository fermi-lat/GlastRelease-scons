
#ifndef digiRootData_CalXtalReadout_H
#define digiRootData_CalXtalReadout_H 1

#include <iostream>
#include <vector>
#include "TObject.h"

#include "CalXtalId.h"

/** @class CalXtalReadout        
* @brief Pulse heights and energy range for both faces of log for Cal.
* Based on CalXtalReadout TDS class in the GlastEvent package written by 
* J. Eric Grove.                                
*            
* @author Heather Kelly
* 
* $Header$
*/
class CalXtalReadout : public TObject  { 
    
public:
    
    CalXtalReadout();

    CalXtalReadout(Char_t rangeP, UShort_t adcP, Char_t rangeM, UShort_t adcM);    
    
    virtual ~CalXtalReadout();

    void initialize(Char_t rangeP, UShort_t adcP, Char_t rangeM, UShort_t adcM);
    
    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    /// retrieve pulse height from specified face
    UShort_t getAdc(CalXtalId::XtalFace face) const;
    
    /// retrieve energy range from specified face
    Char_t getRange(CalXtalId::XtalFace face) const;
    
private:
    
    UShort_t m_adcP;
    UShort_t m_adcM;
    Char_t  m_rangeP;
    Char_t  m_rangeM;
    
    ClassDef(CalXtalReadout,1) // CAL Log Readout
};

#endif
