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
    L1T(UInt_t trigger, UInt_t *digiTriRowBits, UInt_t *trgReqTriRowBits);
    L1T(const L1T& level1);
    virtual ~L1T();

    void initialize(UInt_t trigger, UInt_t *digiTriRowBits, UInt_t *trgReqTriRowBits);
    void setTriggerWordTwo(UInt_t val) { m_triggerWordTwo = val; };
    void setDigiTriRowBits(UInt_t tower, UInt_t bits);
    void setDigiTriRowBits(UInt_t *triRowBits);
    void setTrgReqTriRowBits(UInt_t tower, UInt_t bits);
    void setTrgReqTriRowBits(UInt_t *triRowBits);
    void setPrescale(UInt_t gem, UInt_t glt, UInt_t expired);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    /** @ingroup L1TGroup */
    /*@{*/
    /// Returns full trigger word
    UInt_t getTriggerWord() const { return m_trigger; };
    
    /// Returns the full second word associated with the trigger
    /// Contains the GLT and GEM Engine Numbers
    UInt_t getTriggerWordTwo() const { return m_triggerWordTwo; };

    /// kTRUE indicates that ACD LOW occurred
    //bool getAcdLow() const { return (m_trigger & enums::b_ACDL)!=0;};
    /// kTRUE indicates that ACD HIGH occurred
    bool getAcdHigh() const { return (m_trigger & enums::b_ACDH)!=0;};
    /// kTRUE indicates that 3-in-a-row occurred in the TKR
    bool getTkr3InARow() const { return (m_trigger & enums::b_Track)!=0;};
    /// kTRUE indicates that CAL LOW occurred
    bool getCalLow() const { return (m_trigger & enums::b_LO_CAL)!=0; };
    /// kTRUE indicates that CAL HIGH occurred
    bool getCalHigh() const { return (m_trigger & enums::b_HI_CAL)!=0; };
    /// kTRUE when Ritz Throttle is satisfied
    //bool getThrottle() const { return (m_trigger & enums::b_THROTTLE)!=0; };
    /// kTRUE when ACD ROI is satisfied
    bool getRoi() const { return (m_trigger & enums::b_ROI)!=0; };
	/*@}*/

    /// Corresponds to the GltWord in merit - calculated from digis in TriggerAlg
    UInt_t getGltWord() const { return(m_trigger & 31); };
    /// Corresponds to the GltGemSummary in merit - obtained from GemSummary if real data
    UInt_t getGltGemSummary() const { return ((m_trigger >> enums::GEM_offset) & enums::GEM_mask); };

    /// Corresponds to the GltEngine in merit - the engine number associated with the trigger bits calculated from digis
    Int_t getGltEngine() const { return (m_triggerWordTwo & enums::ENGINE_mask); };
    /// Corresponds to the GltGemEngine in merit - the engine number associated with the GEM summary word
    Int_t getGemEngine() const { return ((m_triggerWordTwo >> enums::ENGINE_offset) & enums::ENGINE_mask); };
    /// Returns kTrue if engine is unset
    bool getEngineUnsetStatus(Int_t engineNum) const { return (engineNum == enums::ENGINE_unset); }

    /// Retrieve the GEM prescale
    UInt_t getGemPrescale() const {return m_gemPrescale; }
 
    /// Retrieve the GLT prescale
    UInt_t getGltPrescale() const {return m_gltPrescale; }

    /// Retrieve prescale expired flag
    UInt_t getPrescaleExpired() const {return m_prescaleExpired; }

    UInt_t getDigiTriRowBits(const Int_t tower) const;
    UInt_t getTrgReqTriRowBits(const Int_t tower) const;

    /// for system tests
    static Int_t number_of_trigger_bits; //!

private:
    /// packed word containing trigger bits
    UInt_t m_trigger;    

    UInt_t m_digiTriRowBits[16];
    UInt_t m_trgReqTriRowBits[16];

    /// packed word containing trigger engine numbers for gem and glt
    /// first 5 bits for glt, mask == 0x1f, second 5 bits for gem, mask == 0x3e00
    UInt_t m_triggerWordTwo;

    /// GEM trigger prescale
    UInt_t  m_gemPrescale;
    /// GLT trigger prescale
    UInt_t  m_gltPrescale;
    /// flag if prescale counter expired on this event
    UInt_t  m_prescaleExpired;

    ClassDef(L1T,8) // Level 1 Trigger information
};

#endif
