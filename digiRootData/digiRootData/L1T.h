// L1T.h
// Dec 05,1999 Richard Dubois Clone from CalHit

#ifndef L1T_H
#define L1T_H

#include "TObject.h"

class L1T : public TObject {
private:
    Int_t    m_trigCount;       // access the 32-bit L1 trigger count    
    Int_t    m_trigTimeValue;   // access the 32-bit the trigger time value (32 bits)    
    Int_t    m_xCapture;        // access the x-trigger REQ capture (32 bits)
    Int_t    m_yCapture;        // access the y-trigger capture (32 bits)    
    UShort_t m_vetoCapture;     // access the veto REQ capture    
    UChar_t  m_deadTimeReason;  // access the dead-time reason (6 bits)    
    UShort_t m_deadTime;        // access the dead-time count (16-bits)
public:

    L1T();
	L1T(Int_t trigCount, Int_t trigTimeValue, Int_t xCapture, Int_t yCapture,
		 UShort_t vetoCapture, UChar_t deadTimeReason, UShort_t deadTime);
    ~L1T();
    Int_t TrigCount(){ return m_trigCount; };
    Int_t TrigTimeValue(){ return m_trigTimeValue; };
    Int_t XCapture(){ return m_xCapture; };
    Int_t YCapture(){ return m_yCapture; };
    UShort_t VetoCapture(){ return m_vetoCapture; };
    UChar_t DeadTimeReason(){ return m_deadTimeReason; };
    UShort_t DeadTime(){ return m_deadTime; };

    void setTrigCount(Int_t newVal) { m_trigCount = newVal; };
    void setTrigTimeValue(Int_t newVal) { m_trigTimeValue = newVal; };
    void setXCapture(Int_t newVal) { m_xCapture = newVal; };
    void setYCapture(Int_t newVal) { m_yCapture = newVal; };
    void setVetoCapture(UShort_t newVal) { m_vetoCapture = newVal; };
    void setDeadTimeReason(UChar_t newVal) { m_deadTimeReason = newVal; };
    void setDeadTime(UShort_t newVal) { m_deadTime = newVal; };

    ClassDef(L1T,2)           // Level 1 Trigger information
};

#endif
