
#ifndef CalLog_H
#define CalLog_H

#include "LogId.h"
#include "TObject.h"

/*! \class CalLog
\brief
The CalLog class contains the information about a single calorimeter  
log.  This class can handle either ADC values for all 4 digitizations, or
a readout that has one "best" range. 
This class implements a LogId member which provides access to information
on the position of the log within the calorimeter.

  Jun 2001 Heather Kelly - renamed from CalHit to CalLog
  Jun 2001 Heather Kelly - remove pointer to LogId
  Jan 2000 Daniel Flath - ROOT HTML comments added
  Dec 1999 Daniel Flath - Rewritten for GLAST
  Oct 25,1999 Richard Dubois Clone from LCD version
*/

class CalLog: public TObject {
public:
    
    CalLog();
    CalLog(LogId &id);
    virtual ~CalLog();
    /// Root >= 3.0 is now const correct for the Compare function
    Int_t Compare(const TObject *obj) const;
    Bool_t IsSortable() const; 
    
    // Enumeration of ADC range values
    // (Should be identical to that found in CalorimeterRawStore.h)
    typedef enum {
        LowEnX4 = 0,
            LowEnX1,
            HiEnX8,
            HiEnX1
    } AdcRange;
    
    typedef enum {
        POS = 0,
            NEG
    } LogFace;
    
    typedef enum
    {
        BESTRANGE = 0,
            ALLRANGE = 2
    } CalTrigMode;
    
private:
/*
ADC word layout
(defined in following enum):
_____________________________________________________
|15|  |  |	|  |  |  |  |  |  |  |04|03|02|01|00|
|__|__|__|__|___|__|__|__|__|__|__|__|__|__|__|__|
|		ADC Value 	        |  ID |PN|RS|
|_______________________________________|_____|__|__|

  ***** NOTES: *****
  
    PN == Pin
    RS == Range Scale
    
      PN & RS = Range
    */
    enum {
        ADC_K_RS = 1,
            ADC_V_RS = 0,
            ADC_M_RS = ((1 << ADC_K_RS) - 1),
            
            ADC_K_PIN = 1,
            ADC_V_PIN = (ADC_K_RS + ADC_V_RS),
            ADC_M_PIN = ((1 << ADC_K_PIN) - 1),
            
            ADC_K_ID = 2,
            ADC_V_ID = (ADC_K_PIN + ADC_V_PIN),
            ADC_M_ID = ((1 << ADC_K_ID) - 1),
            
            ADC_K_VAL = 12,
            ADC_V_VAL = (ADC_K_ID + ADC_V_ID),
            ADC_M_VAL = ((1 << ADC_K_VAL) - 1),
    };
    
    
    
    LogId m_log;       // Log identity info class
    
    //    UShort_t m_AdcValues[2][4];
    UShort_t m_AdcValues[8];
    
    /// Cal readout mode is based on trigger type
    static CalTrigMode m_mode;
    
    
    inline static UShort_t getAdc(UShort_t adcWord) {
        return (adcWord >> ADC_V_VAL) & ADC_M_VAL;
    };
    inline static UShort_t getAdcId(UShort_t adcWord) { 
        return (adcWord >> ADC_V_ID) & ADC_M_ID;
    };
    inline static UShort_t getAdcPinId(UShort_t adcWord) {
        return (adcWord >> ADC_V_PIN) & ADC_M_PIN;
    };
    inline static UShort_t getAdcRangeScale(UShort_t adcWord) {
        return (adcWord >> ADC_V_RS) & ADC_M_RS;
    };
    
    inline static UShort_t getAdcRange(UShort_t adcWord) {
        UShort_t rangeScale = getAdcRangeScale(adcWord);
        UShort_t pinId = getAdcPinId(adcWord);
        return( (pinId << 1 | rangeScale) );
    }
    
    
public:
    
    UShort_t getAdc(LogFace face, short index) const {
        return ( ((face*4+index) < 8) ? getAdc(m_AdcValues[face * 4 + index]) : 0);
    };
    /// Get ADC value for selected face and energy range
    UShort_t getAdcSelectedRange(LogFace face, AdcRange range = LowEnX4) const {
        return ( ((face*4 + range) < 8) ? getAdc(m_AdcValues[face * 4 + range]) : 0);
    };
    //   UShort_t getAdcId(LogFace face, AdcGain range = LowEnX4) const {
    //	return getAdcId(m_AdcValues[end * 4 + gain]);
    //	  return getADCID(m_AdcValues[end][gain]);
    //  };
    //  UShort_t getADCPinId(LogEnd end, ADCGain gain = LowEnX4) const {
    //    return getADCPinId(ADCValues[end * 4 + gain]);
    //	return getADCPinID(m_AdcValues[end][gain]);
    //};
    UShort_t getRange(LogFace face, short index) const {
        return getAdcRange(m_AdcValues[face * 4 + index]);
    };
    
    /// Retrieve readout mode
    inline const CalTrigMode getMode() const { return m_mode; };
    
    
    Bool_t setAdcValue(UShort_t newVal, LogFace face, AdcRange range = LowEnX4);
    Bool_t setAdcId(UShort_t newVal, LogFace face, AdcRange range = LowEnX4);
    Bool_t setAdcPinId(UShort_t newVal, LogFace face, AdcRange range = LowEnX4);
    Bool_t setAdcRangeScale(UShort_t newVal, LogFace face, AdcRange range = LowEnX4);
    inline void setAdcWord(UShort_t newVal, LogFace face, AdcRange range = LowEnX4) {
        m_AdcValues[face * 4 + range] = newVal;
    };
    inline void setTriggerMode(CalTrigMode val) { m_mode = val; };
    
    
    inline LogId* getLogId() { return &m_log; };
    
    
    ClassDef(CalLog,3)	    // Digitization for a single Cal Log
};

#endif

