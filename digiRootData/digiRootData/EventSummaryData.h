#ifndef ROOT_EventSummaryData_H
#define ROOT_EventSummaryData_H

#include "TObject.h"

/** @class EventSummaryData
 * @brief This is the event summary class.
 *
 * The EventSummaryData branch contains:
 * $Header$
 */

class EventSummaryData: public TObject {


public:

    typedef enum {
        GOOD = 0,
        EVTSEQ = 1
    } EventFlags;

    typedef enum {
        GEM = 0,
        OSW = 1,
        ERR = 2,
        DIAG = 3,
        AEM = 4
    }Contribution;


    EventSummaryData();
    EventSummaryData(UInt_t summary) { m_summary = summary; m_flags = 0; };
    EventSummaryData(const EventSummaryData &summary):TObject(summary) {
        m_summary = summary.m_summary;
        m_flags = summary.m_flags;
    }

    virtual ~EventSummaryData();

    void initialize(UInt_t summary);

    void initEventFlags(UInt_t flags) { m_flags = flags; };

    void initTemContribLen(unsigned long *len) {
        unsigned int i;
        for (i = 0; i < 16; i++) { m_temLen[i] = len[i]; }
    }
    void initContribLen(unsigned long *tem, unsigned long gemLen, unsigned long oswLen, unsigned long errLen, unsigned long diagLen, unsigned long aemLen) {

        m_otherContribLen[GEM] = gemLen;
        m_otherContribLen[OSW] = oswLen;
        m_otherContribLen[ERR] = errLen;
        m_otherContribLen[DIAG] = diagLen;
        m_otherContribLen[AEM] = aemLen;
        initTemContribLen(tem);
    }


    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;
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


    /// Some methods we have added
    UInt_t eventSequence() const;
    UInt_t eventFlags() const { return m_flags; };
    Bool_t goodEvent() const { return (m_flags == 0); };
    Bool_t badEvent() const { return (m_flags != 0); };
    Bool_t badEventSequence() const { return (m_flags & EVTSEQ); };
    
    UInt_t temLength(unsigned int tem) { m_temLen[tem]; }
    UInt_t gemLength() const { return m_otherContribLen[GEM]; }
    UInt_t oswLength() const { return m_otherContribLen[OSW]; }
    UInt_t aemLength() const { return m_otherContribLen[AEM]; }
    UInt_t errLength() const { return m_otherContribLen[ERR]; }
    UInt_t diagLength() const { return m_otherContribLen[DIAG]; }

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





