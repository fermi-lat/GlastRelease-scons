
#ifndef digiRootData_CalDigi_H
#define digiRootData_CalDigi_H 1


#include "TObject.h"
#include "CalLogId.h"
#include "CalLogReadout.h"

#ifdef WIN32
using namespace std;
#endif

#include <vector>

/*!
//------------------------------------------------------------------------------
//
// \class   CalDigi        
//  
// \brief Digitizations for Cal                                
//              
// Author:  J. Eric Grove, 23 Feb 2001
//
//------------------------------------------------------------------------------
*/


class CalDigi : public TObject  { 
    
public:
    /// each log end can report four energy ranges
    typedef enum
    {
        LEX8 = 0,
            LEX1,
            HEX8,
            HEX1
    } AdcRange;
    
    /// readout can be either best-of-four energy ranges or all energy ranges
    typedef enum
    {
        BESTRANGE = 0,
            ALLRANGE = 2
    } CalTrigMode;
    
    /// shifts and masks for packed readout of energy range and Adc value
    enum {POS_OFFSET = 14,						// shift for POSitive face
        RANGE_OFFSET = 12, RANGE_MASK = 0x3000,		// energy range bits
        ADC_VAL_MASK = 0xfff};						// Adc value bits
    
    CalDigi();
    
    
    /// Destructor
    virtual ~CalDigi();
    
    
    /// Retrieve readout mode
    const CalTrigMode getMode() const;
    void setMode(CalTrigMode m);
    
    /// Retrieve log identifier
    const CalLogId getPackedId() const;
    void setPackedId(CalLogId id);
    
    void addReadout(CalLogReadout r);
    
    /// Retrieve energy range for selected face and readout
    Char_t getRange(short readoutIndex, CalLogReadout::LogFace face) const;
    
    /// Retrieve pulse height for selected face and readout
    Short_t getAdc(Short_t readoutIndex, CalLogReadout::LogFace face) const;
    
    /// Retrieve ranges and pulse heights from both ends of selected readout
    const CalLogReadout* getLogReadout(Short_t readoutIndex);
    
    /// Retrieve pulse height from selected range
    Short_t getAdcSelectedRange(Char_t range, CalLogReadout::LogFace face) const;
    
private:
    
    /// Cal readout mode is based on trigger type
    CalTrigMode m_mode;
    /// Cal ID
    CalLogId m_logId;
    /// ranges and pulse heights
    vector<CalLogReadout> m_readout;
    
    ClassDef(CalDigi,1)
};




#endif
