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

    EventSummaryData();
    EventSummaryData(UInt_t summary) { m_summary = summary; };
    virtual ~EventSummaryData();

    void initialize(UInt_t summary);

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


private:
    UInt_t m_summary;

    ClassDef(EventSummaryData,1) // Storage for Event Summary Data 
}; 

#endif





