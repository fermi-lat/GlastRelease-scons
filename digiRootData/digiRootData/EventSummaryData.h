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

    EventSummaryData();
    EventSummaryData(UInt_t summary) { m_summary = summary; m_flags = 0; };
    EventSummaryData(const EventSummaryData &summary):TObject(summary) {
        m_summary = summary.m_summary;
        m_flags = summary.m_flags;
    }

    virtual ~EventSummaryData();

    void initialize(UInt_t summary);

    void initEventFlags(UInt_t flags) { m_flags = flags; };

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
    

private:
    UInt_t m_summary;

    UInt_t m_flags;

    ClassDef(EventSummaryData,2) // Storage for Event Summary Data 
}; 

#endif





