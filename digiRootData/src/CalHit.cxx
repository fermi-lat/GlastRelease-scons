///////////////////////////////////////////////////////////////////////////
//                                                                       
// The CalHit class contains the information about a single calorimeter  
// log.  This includes the ADC values for all 4 digitizations, and
// implements a LogID member which provides access to information
// on the position of the log within the calorimeter.
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/CalHit.h"
//#include "TClass.h"

ClassImp(CalHit)

///________________________________________________________________________
CalHit::CalHit() {
    // Default constructor
    m_log = 0;
}
///________________________________________________________________________
CalHit::CalHit(LogID* log) {
    // Create a CalHit object with LogID of log
    m_log = log;
}
//_________________________________________________________________________
CalHit::~CalHit() {
    // Destructor
  //  Clean();
    if (m_log) {
        delete m_log;
        m_log = 0;
    }
}
//_________________________________________________________________________
Int_t CalHit::Compare(const TObject *obj) const {
    if (this == obj) return 0;
    if (CalHit::Class() != obj->IsA()) return -1;
    UInt_t id_this = m_log->getID();
    UInt_t id_hit = ((CalHit*)obj)->getLogID()->getID();
    if (id_this == id_hit)
        return 0;
    else
        return (id_this > id_hit) ? 1 : -1;
}
//_________________________________________________________________________
Bool_t CalHit::IsSortable() const {
  return kTRUE; 
}
//_________________________________________________________________________
//void CalHit::Clean() {
    // Free up memory reserved by member pointers
//    if (m_log)
//        delete m_log;
//}
//_________________________________________________________________________
Bool_t CalHit::setADCValue(UShort_t newVal, LogEnd end, ADCGain gain) {
    // Allows user to set the ADC Value for a particular log end and 
    // digitization.  Returns kTRUE if successful, kFALSE if not.
    if (newVal & ~ADC_M_VAL) {
        return kFALSE;
    } 
    else {
        ADCValues[end * 4 + gain] &= ~(ADC_M_VAL << ADC_V_VAL);
        ADCValues[end * 4 + gain] |= newVal << ADC_V_VAL;
//        ADCValues[end][gain] &= ~(ADC_M_VAL << ADC_V_VAL);
//        ADCValues[end][gain] |= newVal << ADC_V_VAL;
        return kTRUE;
    }
}
//_________________________________________________________________________
Bool_t CalHit::setADCID(UShort_t newVal, LogEnd end, ADCGain gain) {
    // Allows user to set the ADC ID for a particular log end and 
    // digitization.  Returns kTRUE if successful, kFALSE if not.
    if (newVal & ~ADC_M_VAL) {
        return kFALSE;
    } 
    else {
        ADCValues[end * 4 + gain] &= ~(ADC_M_ID << ADC_V_ID);
        ADCValues[end * 4 + gain] |= newVal << ADC_V_ID;
//        ADCValues[end][gain] &= ~(ADC_M_ID << ADC_V_ID);
//        ADCValues[end][gain] |= newVal << ADC_V_ID;
        return kTRUE;
    }
}
//_________________________________________________________________________
Bool_t CalHit::setADCPinID(UShort_t newVal, LogEnd end, ADCGain gain) {
    // Allows user to set the ADC Pin ID for a particular log end and 
    // digitization.  Returns kTRUE if successful, kFALSE if not.
    if (newVal & ~ADC_M_VAL) {
        return kFALSE;
    } 
    else {
        ADCValues[end * 4 + gain] &= ~(ADC_M_PIN << ADC_V_PIN);
        ADCValues[end * 4 + gain] |= newVal << ADC_V_PIN;
//        ADCValues[end][gain] &= ~(ADC_M_PIN << ADC_V_PIN);
//        ADCValues[end][gain] |= newVal << ADC_V_PIN;
        return kTRUE;
    }
}
//_________________________________________________________________________
Bool_t CalHit::setADCRangeScale(UShort_t newVal, LogEnd end, ADCGain gain) {
    // Allows user to set the ADC Range Scale for a particular log end and 
    // digitization.  Returns kTRUE if successful, kFALSE if not.
    if (newVal & ~ADC_M_VAL) {
        return kFALSE;
    } 
    else {
        ADCValues[end * 4 + gain] &= ~(ADC_M_RS << ADC_V_RS);
        ADCValues[end * 4 + gain] |= newVal << ADC_V_RS;
//        ADCValues[end][gain] &= ~(ADC_M_RS << ADC_V_RS);
//        ADCValues[end][gain] |= newVal << ADC_V_RS;
        return kTRUE;
    }
}
