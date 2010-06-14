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
        virtual ~EventSummaryData();
        virtual std::ostream& fillStream(std::ostream &s) const;
        friend std::ostream& operator << (std::ostream &s, const EventSummaryData &obj);

        EventSummaryData() {m_summary=0;};
        EventSummaryData(unsigned int summary) {m_summary=summary; };
        void initialize(unsigned int summary){m_summary=summary;};
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
    };

    inline EventSummaryData::~EventSummaryData(){
    }


    inline std::ostream& EventSummaryData::fillStream(std::ostream &s) const{
        s << "EventSummary\n"
          << "calStrobe: " << calStrobe() << " tag: " << tag()
         << "\nTACK: " << TACK() << " readout4: " << readout4()
          << "\nzeroSuppress: " << zeroSuppress() << " marker: " << marker()
          << "\nerror: " << error() << " diagnostic: " << diagnostic()
          << "\neventNumber: " << eventNumber() << " trgParityError: " << trgParityError() << "\n";
        return s;
    }


    inline std::ostream& operator<<(std::ostream &s, const EventSummaryData &obj){
        return obj.fillStream(s);
    }

}//namespace LdfEvent 


#endif
