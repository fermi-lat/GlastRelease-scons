
#ifndef digiRootData_CalDigi_H
#define digiRootData_CalDigi_H 1


#include "TObject.h"
#include "CalXtalId.h"
#include "CalXtalReadout.h"

#ifdef WIN32
using namespace std;
#endif

#include <vector>

/** @class CalDigi
 * @brief Digitizations for CAL                                
 *
 * @author J. Eric Grove
 * $Header$
*/

class CalDigi : public TObject  { 
    
public:
    
    /// shifts and masks for packed readout of energy range and Adc value
    enum {POS_OFFSET = 14,						// shift for POSitive face
        RANGE_OFFSET = 12, RANGE_MASK = 0x3000, // energy range bits
        ADC_VAL_MASK = 0xfff};						// Adc value bits
    
    CalDigi();
    
    virtual ~CalDigi();
    
    /// Retrieve readout mode
    const CalXtalId::CalTrigMode getMode() const;
    void setMode(CalXtalId::CalTrigMode m);
    
    /// Retrieve log identifier
    const CalXtalId getPackedId() const;
    void setPackedId(CalXtalId id);
    
    void addReadout(CalXtalReadout r);
    
    /// Retrieve energy range for selected face and readout
    Char_t getRange(short readoutIndex, CalXtalId::XtalFace face) const;
    
    /// Retrieve pulse height for selected face and readout
    Short_t getAdc(Short_t readoutIndex, CalXtalId::XtalFace face) const;
    
    /// Retrieve ranges and pulse heights from both ends of selected readout
    const CalXtalReadout* getLogReadout(Short_t readoutIndex);
    
    /// Retrieve pulse height from selected range
    Short_t getAdcSelectedRange(Char_t range, CalXtalId::XtalFace face) const;
    
private:
    
    /// Cal readout mode is based on trigger type
    CalXtalId::CalTrigMode m_mode;
    /// Cal Xtal Id
    CalXtalId m_xtalId;
    /// ranges and pulse heights
    vector<CalXtalReadout> m_readout;
    
    ClassDef(CalDigi,1) // CAL Digitization Class
};

#endif
