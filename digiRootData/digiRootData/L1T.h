#ifndef ROOT_L1T_H
#define ROOT_L1T_H

#include "TObject.h"

#include "enums/TriggerBits.h"

/** @class L1T
 * @brief Store the Level One Trigger information
 *
 * $Header$
*/

class L1T : public TObject {
public:
    
    L1T();
    L1T(UInt_t trigger, UInt_t *triRowBits);
    L1T(const L1T& level1);
    virtual ~L1T();

    void initialize(UInt_t trigger, UInt_t *triRowBits);
    void setTriRowBits(UInt_t tower, UInt_t bits);
    void setTriRowBits(UInt_t *triRowBits);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    /** @ingroup L1TGroup */
    /*@{*/
    UInt_t getTriggerWord() const { return m_trigger; };

    /// kTRUE indicates that ACD LOW occurred
    bool getAcdLow() const { return m_trigger & enums::b_ACDL;};
    /// kTRUE indicates that ACD HIGH occurred
    bool getAcdHigh() const { return m_trigger & enums::b_ACDH;};
    /// kTRUE indicates that 3-in-a-row occurred in the TKR
    bool getTkr3InARow() const { return m_trigger & enums::b_Track;};
    /// kTRUE indicates that CAL LOW occurred
    bool getCalLow() const { return m_trigger & enums::b_LO_CAL; };
    /// kTRUE indicates that CAL HIGH occurred
    bool getCalHigh() const { return m_trigger & enums::b_HI_CAL; };
    /// kTRUE when Ritz Throttle is satisfied
    bool getThrottle() const { return m_trigger & enums::b_THROTTLE; };
	/*@}*/

    UInt_t getTriRowBits(const Int_t tower) const;

    /// for system tests
    static Int_t number_of_trigger_bits; //!

private:
    /// packed word containing trigger bits
    UInt_t    m_trigger;    

    UInt_t m_triRowBits[16];


    ClassDef(L1T,5) // Level 1 Trigger information
};

#endif
