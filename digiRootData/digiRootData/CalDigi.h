
#ifndef digiRootData_CalDigi_H
#define digiRootData_CalDigi_H 1

#include "TObject.h"
#include "TClonesArray.h"
#include "commonRootData/idents/CalXtalId.h"
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
    
    /** @defgroup CalDigiGroup CalDigi End-User Interface */
    /*@{*/
    /// shifts and masks for packed readout of energy range and Adc value
    enum {
        POS_OFFSET = 14,     // shift for POSitive face
        RANGE_OFFSET = 12, 
        RANGE_MASK = 0x3000, // energy range bits
        ADC_VAL_MASK = 0xfff
    };
    /*@}*/
    
    CalDigi();
    
    virtual ~CalDigi();

    void initialize(CalXtalId::CalTrigMode m, const CalXtalId &id);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    /** @ingroup CalDigiGroup */
    /*@{*/
    /// Retrieve readout mode - BEST or ALLRANGE
    const CalXtalId::CalTrigMode getMode() const { return m_mode; };
    
    /// Retrieve log identifier
    const CalXtalId getPackedId() const { return m_xtalId; };
    
    /// Returns the whole collection of readouts
    const CalXtalReadout* getReadoutCol() const;

    /// Number of readouts in this CalDigi (1 or 4)
    const UInt_t getNumReadouts() const { 
        if (m_readoutCol) 
            return m_readoutCol->GetEntries(); 
        else 
            return m_numReadouts;
    }
    
    /// Retrieve energy range for selected face and readout
    Char_t getRange(UShort_t readoutIndex, CalXtalId::XtalFace face) const;
    
    /// Retrieve pulse height for selected face and readout
    Short_t getAdc(UShort_t readoutIndex, CalXtalId::XtalFace face) const;
    
    /// Retrieve ranges and pulse heights from both ends of selected readout
    const CalXtalReadout* getXtalReadout(UShort_t readoutIndex) const;
    
    /// Retrieve pulse height from selected range
    Short_t getAdcSelectedRange(Char_t range, CalXtalId::XtalFace face) const;
    /*@}*/

    const CalXtalReadout* addReadout(Char_t rangeP, UShort_t adcP, Char_t rangeM, UShort_t adcM);

private:
    
    /// Cal readout mode is based on trigger type
    CalXtalId::CalTrigMode m_mode;
    /// number of readouts actually filled in m_readoutArr
    UInt_t m_numReadouts;
    /// Cal Xtal Id
    CalXtalId m_xtalId;
    /// collection of ranges and pulse heights
    // Allow for backward compatibilty by retaining the old TClonesArray*
    TClonesArray *m_readoutCol;
    CalXtalReadout m_readoutArr[4];
    
    ClassDef(CalDigi,2) // CAL Digitization Class
};

#endif
