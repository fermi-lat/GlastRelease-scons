#ifndef LDFCONVERTER_EventSummaryDATA_H
#define LDFCONVERTER_EventSummaryDATA_H

#include <iostream>
#include <vector>

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"

#include "Event/TopLevel/Definitions.h"
#include "Event/TopLevel/EventModel.h"

#include "eventSummary.h" // from EBF library

/**
* @class EventSummaryData
* @brief TDS for storing the event summary data 
*/
extern const CLID& CLID_LdfEventSummaryData;
namespace LdfEvent{

    class EventSummaryData : public DataObject {
    public:
        /// EVTSEQ set means that the event sequence number
        /// constructed from the eventNumber and tag in the EventSummary
        /// are either inconsistent across contributions in the evnet
        /// or the event sequence across events is not monotonically increasing
        typedef enum {
            GOOD = 0,
            EVTSEQ = 1
        } EventFlags ;

        virtual ~EventSummaryData();
        virtual std::ostream& fillStream(std::ostream &s) const;
        friend std::ostream& operator << (std::ostream &s, const EventSummaryData &obj);

        EventSummaryData() {m_summary=0; m_flags = 0;};
        EventSummaryData(unsigned int summary) {m_summary=summary; m_flags = 0;};
        void initialize(unsigned int summary){m_summary=summary;};
        void initEventFlags(unsigned int flags) { m_flags = flags; };
        
        unsigned int eventFlags() const { return m_flags; };
        bool goodEvent() const { return (m_flags == 0); };
        bool badEventSeq() const { return (m_flags && EVTSEQ); };
        bool badEvent() const { return (m_flags != 0); } ;
        unsigned long eventSequence() const {
            unsigned eventNumber = EventSummary::eventNumber(m_summary);
            unsigned tag = EventSummary::tag(m_summary);
            return ((eventNumber << 2) | tag);
        }
 
        unsigned int summary() { return m_summary; };

        unsigned int calStrobe() const {return EventSummary::calStrobe(m_summary);};
        unsigned int tag() const {return EventSummary::tag(m_summary);};
        unsigned int TACK() const {return EventSummary::TACK(m_summary);};
        unsigned int readout4() const {return EventSummary::readout4(m_summary);};
        unsigned int zeroSuppress() const {return EventSummary::zeroSuppress(m_summary);};
        unsigned int marker() const {return EventSummary::marker(m_summary);};
        unsigned int error() const {return EventSummary::error(m_summary);};
        unsigned int diagnostic() const {return EventSummary::diagnostic(m_summary);};
        unsigned int eventNumber() const {return EventSummary::eventNumber(m_summary);};
        unsigned int trgParityError() const {return EventSummary::trgParityError(m_summary);}; 

    private:
    unsigned int m_summary;
    /// Event Flags - used to denote bad events
    /// flags = 0 is "good"
    unsigned int m_flags;
    };

    inline EventSummaryData::~EventSummaryData(){
    }


    inline std::ostream& EventSummaryData::fillStream(std::ostream &s) const{
        s << "EventSummary\n"
          << "calStrobe: " << calStrobe() << " tag: " << tag()
         << "\nTACK: " << TACK() << " readout4: " << readout4()
          << "\nzeroSuppress: " << zeroSuppress() << " marker: " << marker()
          << "\nerror: " << error() << " diagnostic: " << diagnostic()
          << "\neventNumber: " << eventNumber() << " trgParityError: " << trgParityError() << "\n"
          << "Event Flags: " << m_flags << "\n";
        return s;
    }


    inline std::ostream& operator<<(std::ostream &s, const EventSummaryData &obj){
        return obj.fillStream(s);
    }

}//namespace LdfEvent 


#endif
