#ifndef ROOT_EventSummaryData_H
#define ROOT_EventSummaryData_H

#include "TObject.h"
#include "enums/EventFlags.h"

/** @class EventSummaryData
 * @brief This is the event summary class.
 *
 * The EventSummaryData branch contains:
 * $Header$
 */

class EventSummaryData: public TObject {


public:

    /**@defgroup EventSummaryDataGroup EventSummaryData End-User Interface */
	/*@{*/
    typedef enum {
        GEM = 0,
        OSW = 1,
        AEM = 2
    }Contribution;
	/*@}*/

    EventSummaryData();
    EventSummaryData(UInt_t summary) { Clear(); m_summary = summary; m_flags = 0; };
    EventSummaryData(const EventSummaryData &summary):TObject(summary) {
        m_summary = summary.m_summary;
        m_flags = summary.m_flags;
        m_otherContribLen[GEM] = summary.m_otherContribLen[GEM];
        m_otherContribLen[OSW] = summary.m_otherContribLen[OSW];
        m_otherContribLen[AEM] = summary.m_otherContribLen[AEM];
        int i;
        for (i = 0; i < 16; i++) { 
           m_temLen[i] = summary.m_temLen[i]; 
           m_errLen[i] = summary.m_errLen[i];
           m_diagLen[i] = summary.m_diagLen[i];
        }
        m_eventSizeInBytes = summary.m_eventSizeInBytes;
    }

    virtual ~EventSummaryData();

    void initialize(UInt_t summary);

    //void initEventSequence(UInt_t evtSeq) { m_eventSequence = evtSeq; };

    void initEventFlags(UInt_t flags) { m_flags = flags; };

    void initEventSizeInBytes(ULong_t size) { m_eventSizeInBytes = size; }

    void initTemContribLen(unsigned int *len, unsigned int *diag, unsigned int *err) {
        unsigned int i;
        for (i = 0; i < 16; i++) { 
            m_temLen[i] = len[i]; 
            m_errLen[i] = err[i];
            m_diagLen[i] = diag[i];
        }
    }
    void initContribLen(unsigned int *tem, unsigned int gemLen, unsigned int oswLen, unsigned int *errLen, unsigned int *diagLen, unsigned int aemLen) {

        m_otherContribLen[GEM] = gemLen;
        m_otherContribLen[OSW] = oswLen;
        m_otherContribLen[AEM] = aemLen;
        initTemContribLen(tem, diagLen, errLen);
    }


    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;
    
	/** @ingroup EventSummaryDataGroup */
	/*@{*/
	UInt_t summary() { return m_summary; };

    UInt_t calStrobe();
    UInt_t tag();
    UInt_t TACK();
    UInt_t readout4();
    UInt_t zeroSuppress();
    UInt_t marker();
    UInt_t error();
    UInt_t diagnostic();
    UInt_t eventNumber();
    // renamed to promote the other trgParityError method
    UInt_t getTrgParityError();


    ULong_t eventSizeInBytes() const { return m_eventSizeInBytes; };

    UInt_t eventSequence() const { return m_eventSequence; };
    UInt_t eventFlags() const { return m_flags; };
    /// Returns true if the flags member is zero
    Bool_t goodEvent() const { return (m_flags == 0); };
    /// Returns true if the flag member is non-zero
    Bool_t badEvent() const { return (m_flags != 0); };
    /// Checks of event sequence bit is set in the event flags
    Bool_t badEventSequence() const { return (m_flags & enums::EVTSEQ)!=0; };
    /// Checks the TkrRecon bit, if set an error occurred during TkrRecon
    Bool_t badTkrRecon() const { return (m_flags & enums::TKRRECON)!=0; };
    Bool_t packetError() const { return (m_flags & enums::PACKETERROR)!=0; };
    Bool_t temError() const { return (m_flags & enums::SUMMARYERROR)!=0; };
    // old name for temError method
    Bool_t errorEventSummary() const { return (m_flags & enums::SUMMARYERROR)!=0; };
    Bool_t trgParityError() const { return (m_flags & enums::TRGPARITYERROR)!=0; };
    
    /// Returns the length in bytes of the TEM contribution identified by a value in [0,15]
    UInt_t temLength(unsigned int tem) const { return m_temLen[tem]; }
    const UInt_t* temLength() const { return m_temLen; };
    /// Returns the length in bytes of the GEM contribution
    UInt_t gemLength() const { return m_otherContribLen[GEM]; }
    /// Returns the length in bytes of the OSW contribution
    UInt_t oswLength() const { return m_otherContribLen[OSW]; }
    /// Returns the length in bytes of the AEM contribution
    UInt_t aemLength() const { return m_otherContribLen[AEM]; }
    /// Returns the length in bytes of the ERROR contribution
    UInt_t errLength(unsigned int tem) const { return m_errLen[tem]; }
    const UInt_t* errLength() const { return m_errLen; };
    /// Returns the length in bytes of the DIAGNOSTIC contribution
    UInt_t diagLength(unsigned int tem) const { return m_diagLen[tem]; }
    const UInt_t* diagLength() const { return m_diagLen; };
    /*@}*/

private:
    UInt_t m_summary;

    UInt_t m_flags;

    /// Stores the lengths from each TEM contribution
    UInt_t m_temLen[16];
    UInt_t m_diagLen[16], m_errLen[16];
    /// Stores the lengths of the "other" contributions
    /// GEM, AEM,, OSW
    UInt_t m_otherContribLen[3];

    UInt_t m_eventSequence;
    ULong_t m_eventSizeInBytes;

    ClassDef(EventSummaryData,6) // Storage for Event Summary Data 
}; 

#endif





