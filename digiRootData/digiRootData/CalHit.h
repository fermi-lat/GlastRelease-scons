
#ifndef CALHIT_H
#define CALHIT_H

#include "LogID.h"
#include "TObject.h"

/*! \class CalHit
\brief
 Jun 2001 Heather Kelly - remove pointer to LogID
 Jan 2000 Daniel Flath - ROOT HTML comments added
 Dec 1999 Daniel Flath - Rewritten for GLAST
 Oct 25,1999 Richard Dubois Clone from LCD version
*/

class CalHit: public TObject {
private:
    /*
    ADC word layout
    (defined in following enum):
     _______________________________________________
    |15|  |  |	|  |  |  |  |  |  |  |04|03|02|01|00|
    |__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
    |		  Value 		|  ID |PN|RS|
    |___________________________________|_____|__|__|

    ***** NOTES: *****

    PN == Pin
    RS == Range Scale
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

    LogID m_log;       // Log identity info class

//    UShort_t ADCValues[2][4];
    UShort_t ADCValues[8];

    inline static UShort_t getADCValue(UShort_t adcWord) {
	return (adcWord >> ADC_V_VAL) & ADC_M_VAL;
    };
    inline static UShort_t getADCID(UShort_t adcWord) { 
	return (adcWord >> ADC_V_ID) & ADC_M_ID;
    };
    inline static UShort_t getADCPinID(UShort_t adcWord) {
	return (adcWord >> ADC_V_PIN) & ADC_M_PIN;
    };
    inline static UShort_t getADCRangeScale(UShort_t adcWord) {
	return (adcWord >> ADC_V_RS) & ADC_M_RS;
    };
public:

    // Enumeration of ADC gain values
    // (Should be identical to that found in CalorimeterRawStore.h)
    typedef enum {
	LowEnX4 = 0,
	LowEn,
	HiEnX8,
	HiEn
    } ADCGain;

    typedef enum {
	Pos = 0,
	Neg
    } LogEnd;


    UShort_t getADCValue(LogEnd end, ADCGain gain = LowEnX4) const {
	return getADCValue(ADCValues[end * 4 + gain]);
//	  return getADCValue(ADCValues[end][gain]);
    };
    UShort_t getADCID(LogEnd end, ADCGain gain = LowEnX4) const {
	return getADCID(ADCValues[end * 4 + gain]);
//	  return getADCID(ADCValues[end][gain]);
    };
    UShort_t getADCPinID(LogEnd end, ADCGain gain = LowEnX4) const {
      return getADCPinID(ADCValues[end * 4 + gain]);
//	return getADCPinID(ADCValues[end][gain]);
    };
    UShort_t getADCRangeScale(LogEnd end, ADCGain gain = LowEnX4) const {
      return getADCRangeScale(ADCValues[end * 4 + gain]);
//	return getADCRangeScale(ADCValues[end][gain]);
    };

    Bool_t setADCValue(UShort_t newVal, LogEnd end, ADCGain gain = LowEnX4);
    Bool_t setADCID(UShort_t newVal, LogEnd end, ADCGain gain = LowEnX4);
    Bool_t setADCPinID(UShort_t newVal, LogEnd end, ADCGain gain = LowEnX4);
    Bool_t setADCRangeScale(UShort_t newVal, LogEnd end, ADCGain gain = LowEnX4);
    inline void setADCWord(UShort_t newVal, LogEnd end, ADCGain gain = LowEnX4) {
	ADCValues[end * 4 + gain] = newVal;
    };

    CalHit();
    //CalHit(LogID* log);
    virtual ~CalHit();
    inline LogID* getLogID() { return &m_log; };
   // void Clean();
    
    /// Root >= 3.0 is now const correct for the Compare function
    Int_t Compare(const TObject *obj) const;
    Bool_t IsSortable() const; 

    ClassDef(CalHit,3)	    // Information for a single CsI Log
};

#endif

