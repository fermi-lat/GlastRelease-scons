
#ifndef L1T_H
#define L1T_H

#include "TObject.h"

/*! \class L1T
\brief Store the Level One Trigger information
 
\li Dec 05,1999 Richard Dubois Clone from CalHit
*/

class L1T : public TObject {
private:
    UInt_t    m_trigCount;       // access the 32-bit L1 trigger count    
    UInt_t    m_trigTimeValue;   // access the 32-bit the trigger time value (32 bits)    
    UInt_t    m_xCapture;        // access the x-trigger REQ capture (32 bits)
    UInt_t    m_yCapture;        // access the y-trigger capture (32 bits)    
    UShort_t m_vetoCapture;     // access the veto REQ capture    
    UInt_t m_deadTimeAndCause;  // original 32 bit word that contains both dead time and cause
    UChar_t  m_deadTimeCause;  // access the dead-time reason (6 bits)    
    UShort_t m_deadTime;        // access the dead-time count (16-bits)
public:
    enum {
        DeadTimeMask = 0x00002fff,
        DeadTimeCauseMask = 0x001fc000,
        DeadTimeCauseShift = 14
    };

    L1T();
	L1T(Int_t trigCount, Int_t trigTimeValue, Int_t xCapture, Int_t yCapture,
		 UShort_t vetoCapture, UChar_t deadTimeCause, UShort_t deadTime);
    virtual ~L1T();
    void Clean();
    Int_t getTrigCount(){ return m_trigCount; };
    Int_t getTrigTime(){ return m_trigTimeValue; };
    Int_t getXCapture(){ return m_xCapture; };
    Int_t getYCapture(){ return m_yCapture; };
    UShort_t getVetoCapture(){ return m_vetoCapture; };
    UInt_t getDeadTimeAndCause() { return m_deadTimeAndCause; };
    UChar_t getDeadTimeCause(){ return ((m_deadTimeAndCause & DeadTimeCauseMask) >> DeadTimeCauseShift); };
    UShort_t getDeadTime(){ return (m_deadTimeAndCause & DeadTimeMask); };

    void setTrigCount(Int_t newVal) { m_trigCount = newVal; };
    void setTrigTime(Int_t newVal) { m_trigTimeValue = newVal; };
    void setXCapture(Int_t newVal) { m_xCapture = newVal; };
    void setYCapture(Int_t newVal) { m_yCapture = newVal; };
    void setVetoCapture(UShort_t newVal) { m_vetoCapture = newVal; };
    void setDeadTimeAndCause(UInt_t newVal) { m_deadTimeAndCause = newVal; };
    //void setDeadTime(UShort_t newVal) { m_deadTime = newVal; };

    ClassDef(L1T,3)           // Level 1 Trigger information
};

#endif
