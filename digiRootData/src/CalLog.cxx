
#include "digiRootData/CalLog.h"

ClassImp(CalLog)

CalLog::CalTrigMode CalLog::m_mode;

///________________________________________________________________________
CalLog::CalLog() {
    // Default constructor
    m_mode = ALLRANGE;
    m_log.setTag(UInt_t(0));
}
///________________________________________________________________________
CalLog::~CalLog() {
    // Destructor
}
//_________________________________________________________________________
Int_t CalLog::Compare(const TObject *obj) const {
    if (this == obj) return 0;
    if (CalLog::Class() != obj->IsA()) return -1;
    UInt_t id_this = m_log.getId();
    UInt_t id_hit = ((CalLog*)obj)->getLogId()->getId();
    if (id_this == id_hit)
        return 0;
    else
        return (id_this > id_hit) ? 1 : -1;
}
//_________________________________________________________________________
Bool_t CalLog::IsSortable() const {
    return kTRUE; 
}
//_________________________________________________________________________
Bool_t CalLog::setAdcValue(UShort_t newVal, LogFace face, AdcRange range) {
    // Allows user to set the ADC Value for a particular log end and 
    // digitization.  Returns kTRUE if successful, kFALSE if not.
    if (newVal & ~ADC_M_VAL) {
        return kFALSE;
    } 
    else {
        m_AdcValues[face * 4 + range] &= ~(ADC_M_VAL << ADC_V_VAL);
        m_AdcValues[face * 4 + range] |= newVal << ADC_V_VAL;
        //	  m_AdcValues[end][gain] &= ~(ADC_M_VAL << ADC_V_VAL);
        //	  m_AdcValues[end][gain] |= newVal << ADC_V_VAL;
        return kTRUE;
    }
}
//_________________________________________________________________________
Bool_t CalLog::setAdcId(UShort_t newVal, LogFace face, AdcRange range) {
    // Allows user to set the ADC ID for a particular log end and 
    // digitization.  Returns kTRUE if successful, kFALSE if not.
    if (newVal & ~ADC_M_VAL) {
        return kFALSE;
    } 
    else {
        m_AdcValues[face * 4 + range] &= ~(ADC_M_ID << ADC_V_ID);
        m_AdcValues[face * 4 + range] |= newVal << ADC_V_ID;
        //	  m_AdcValues[end][gain] &= ~(ADC_M_ID << ADC_V_ID);
        //	  m_AdcValues[end][gain] |= newVal << ADC_V_ID;
        return kTRUE;
    }
}
//_________________________________________________________________________
Bool_t CalLog::setAdcPinId(UShort_t newVal, LogFace face, AdcRange range) {
    // Allows user to set the ADC Pin ID for a particular log end and 
    // digitization.  Returns kTRUE if successful, kFALSE if not.
    if (newVal & ~ADC_M_VAL) {
        return kFALSE;
    } 
    else {
        m_AdcValues[face * 4 + range] &= ~(ADC_M_PIN << ADC_V_PIN);
        m_AdcValues[face * 4 + range] |= newVal << ADC_V_PIN;
        //	  m_AdcValues[end][gain] &= ~(ADC_M_PIN << ADC_V_PIN);
        //	  m_AdcValues[end][gain] |= newVal << ADC_V_PIN;
        return kTRUE;
    }
}
//_________________________________________________________________________
Bool_t CalLog::setAdcRangeScale(UShort_t newVal, LogFace face, AdcRange range) {
    // Allows user to set the ADC Range Scale for a particular log end and 
    // digitization.  Returns kTRUE if successful, kFALSE if not.
    if (newVal & ~ADC_M_VAL) {
        return kFALSE;
    } 
    else {
        m_AdcValues[face * 4 + range] &= ~(ADC_M_RS << ADC_V_RS);
        m_AdcValues[face * 4 + range] |= newVal << ADC_V_RS;
        //	  m_AdcValues[end][gain] &= ~(ADC_M_RS << ADC_V_RS);
        //	  m_AdcValues[end][gain] |= newVal << ADC_V_RS;
        return kTRUE;
    }
}
