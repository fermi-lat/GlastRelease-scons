
#ifndef digiRootData_CalDigi_H
#define digiRootData_CalDigi_H 1

#include "TObject.h"
#include "TClonesArray.h"
#include "CalXtalId.h"
#include "CalXtalReadout.h"

/** @class CalDigi
 * @brief Digitizations for CAL
 *
 * Based on the CalDigi TDS class in GlastEvent written by J. Eric Grove.
 *
 * @author Heather Kelly
 * $Header$
*/

class CalDigi : public TObject  { 
    
public:
    
    /// shifts and masks for packed readout of energy range and Adc value
    enum {
        POS_OFFSET = 14,     // shift for POSitive face
        RANGE_OFFSET = 12, 
        RANGE_MASK = 0x3000, // energy range bits
        ADC_VAL_MASK = 0xfff
    };
    
    CalDigi();
    
    virtual ~CalDigi();

    void initialize(CalXtalId::CalTrigMode m, const CalXtalId &id);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    /// Retrieve readout mode
    const CalXtalId::CalTrigMode getMode() const { return m_mode; };
    
    /// Retrieve log identifier
    const CalXtalId getPackedId() const { return m_xtalId; };
    
    const CalXtalReadout* addReadout(Char_t rangeP, UShort_t adcP, Char_t rangeM, UShort_t adcM);
    const TClonesArray* getReadoutCol() const { return m_readoutCol; };
    
    /// Retrieve energy range for selected face and readout
    Char_t getRange(short readoutIndex, CalXtalId::XtalFace face) const;
    
    /// Retrieve pulse height for selected face and readout
    Short_t getAdc(Short_t readoutIndex, CalXtalId::XtalFace face) const;
    
    /// Retrieve ranges and pulse heights from both ends of selected readout
    const CalXtalReadout* getXtalReadout(Short_t readoutIndex);
    
    /// Retrieve pulse height from selected range
    Short_t getAdcSelectedRange(Char_t range, CalXtalId::XtalFace face) const;
    
private:
    
    /// Cal readout mode is based on trigger type
    CalXtalId::CalTrigMode m_mode;
    UInt_t m_numXtals;
    /// Cal Xtal Id
    CalXtalId m_xtalId;
    /// collection of ranges and pulse heights
    TClonesArray *m_readoutCol;
    
    ClassDef(CalDigi,1) // CAL Digitization Class
};

#endif
