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
        ERR = 2,
        DIAG = 3,
        AEM = 4
    }Contribution;
	/*@}*/

    EventSummaryData();
    EventSummaryData(UInt_t summary) { Clear(); m_summary = summary; m_flags = 0; };
    EventSummaryData(const EventSummaryData &summary):TObject(summary) {
        m_summary = summary.m_summary;
        m_flags = summary.m_flags;
        m_otherContribLen[GEM] = summary.m_otherContribLen[GEM];
        m_otherContribLen[OSW] = summary.m_otherContribLen[OSW];
        m_otherContribLen[ERR] = summary.m_otherContribLen[ERR];
        m_otherContribLen[DIAG] = summary.m_otherContribLen[DIAG];
        m_otherContribLen[AEM] = summary.m_otherContribLen[AEM];
        int i;
        for (i = 0; i < 16; i++) { m_temLen[i] = summary.m_temLen[i]; }
    }

    virtual ~EventSummaryData();

    void initialize(UInt_t summary);

    void initEventFlags(UInt_t flags) { m_flags = flags; };

    void initTemContribLen(unsigned int *len) {
        unsigned int i;
        for (i = 0; i < 16; i++) { m_temLen[i] = len[i]; }
    }
    void initContribLen(unsigned int *tem, unsigned int gemLen, unsigned int oswLen, unsigned int errLen, unsigned int diagLen, unsigned int aemLen) {

        m_otherContribLen[GEM] = gemLen;
        m_otherContribLen[OSW] = oswLen;
        m_otherContribLen[ERR] = errLen;
        m_otherContribLen[DIAG] = diagLen;
        m_otherContribLen[AEM] = aemLen;
        initTemContribLen(tem);
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
    UInt_t trgParityError();


    UInt_t eventSequence() const;
    UInt_t eventFlags() const { return m_flags; };
    /// Returns true if the flags member is zero
    Bool_t goodEvent() const { return (m_flags == 0); };
    /// Returns true if the flag member is non-zero
    Bool_t badEvent() const { return (m_flags != 0); };
    /// Checks of event sequence bit is set in the event flags
    Bool_t badEventSequence() const { return (m_flags & enums::EVTSEQ); };
    /// Checks the TkrRecon bit, if set an error occurred during TkrRecon
    Bool_t badTkrRecon() const { return (m_flags & enums::TKRRECON); };
    Bool_t packetError() const { return (m_flags & enums::PACKETERROR); };
    Bool_t errorEventSummary() const { return (m_flags & enums::SUMMARYERROR); };
    
    /// Returns the length in bytes of the TEM contribution identified by a value in [0,15]
    UInt_t temLength(unsigned int tem) { return m_temLen[tem]; }
    /// Returns the length in bytes of the GEM contribution
    UInt_t gemLength() const { return m_otherContribLen[GEM]; }
    /// Returns the length in bytes of the OSW contribution
    UInt_t oswLength() const { return m_otherContribLen[OSW]; }
    /// Returns the length in bytes of the AEM contribution
    UInt_t aemLength() const { return m_otherContribLen[AEM]; }
    /// Returns the length in bytes of the ERROR contribution
    UInt_t errLength() const { return m_otherContribLen[ERR]; }
    /// Returns the length in bytes of the DIAGNOSTIC contribution
    UInt_t diagLength() const { return m_otherContribLen[DIAG]; }
    /*@}*/

private:
    UInt_t m_summary;

    UInt_t m_flags;

    /// Stores the lengths from each TEM contribution
    UInt_t m_temLen[16];
    /// Stores the lengths of the "other" contributions
    /// GEM, AEM, Error, DIAGNOSTIC, OSW
    UInt_t m_otherContribLen[5];

    ClassDef(EventSummaryData,3) // Storage for Event Summary Data 
}; 

#endif





