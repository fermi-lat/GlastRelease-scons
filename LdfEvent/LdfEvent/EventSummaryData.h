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
       /// TKRRECON set means there was some error condition during TkrRecon
        typedef enum {
            GOOD = 0,
            EVTSEQ = 1,
            TKRRECON = 2
        } EventFlags ;


        typedef enum {
            GEM = 0,
            OSW = 1,
            ERR = 2, 
            DIAG = 3,
            AEM = 4
        }Contribution;

        virtual ~EventSummaryData();
        virtual std::ostream& fillStream(std::ostream &s) const;
        friend std::ostream& operator << (std::ostream &s, const EventSummaryData &obj);

        EventSummaryData() { clear(); };
        EventSummaryData(unsigned int summary) {clear(); m_summary=summary; };

        /// Retrieve reference to class definition structure
        virtual const CLID& clID() const   { return EventSummaryData::classID(); }
        static const CLID& classID()       { return CLID_LdfEventSummaryData; }


        void clear() { 
            m_summary = 0; m_flags = 0;
            unsigned int i;
            for (i = 0; i<16; i++) m_temLen[i] = 0;
            m_otherContribLen[GEM] = 0;
            m_otherContribLen[AEM] = 0;
            m_otherContribLen[OSW] = 0;
            m_otherContribLen[ERR] = 0;
            m_otherContribLen[DIAG] = 0;
        }

        void initialize(unsigned int summary){m_summary=summary;};
        void initEventFlags(unsigned int flags) { m_flags = flags; };
        void initTemContribLen(unsigned int *len) {
            unsigned int i;
            for (i = 0; i < 16; i++) { m_temLen[i] = len[i]; }
        }
        void initContribLen(unsigned int *tem, unsigned int gemLen, unsigned int oswLen,
            unsigned int errLen, unsigned int diagLen, unsigned int aemLen) {

            m_otherContribLen[GEM] = gemLen;
            m_otherContribLen[OSW] = oswLen;
            m_otherContribLen[ERR] = errLen;
            m_otherContribLen[DIAG] = diagLen;
            m_otherContribLen[AEM] = aemLen;
            initTemContribLen(tem);
        }
        
        unsigned int temLength(unsigned int tem) { return m_temLen[tem]; }
        unsigned int gemLength() const { return m_otherContribLen[GEM]; }
        unsigned int oswLength() const { return m_otherContribLen[OSW]; }
        unsigned int aemLength() const { return m_otherContribLen[AEM]; }
        unsigned int errLength() const { return m_otherContribLen[ERR]; }
        unsigned int diagLength() const { return m_otherContribLen[DIAG]; }

        void setEvtSeqBit()  { m_flags |= EVTSEQ; };
        void setTkrReconBit() { m_flags |= TKRRECON; };
        unsigned int eventFlags() const { return m_flags; };
        bool goodEvent() const { return (m_flags == 0); };
        bool badEventSeq() const { return (m_flags && EVTSEQ); };
        bool badTkrRecon() const { return (m_flags && TKRRECON); };
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
    unsigned int m_temLen[16];
    unsigned int m_otherContribLen[5];
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
          << "Event Flags: " << m_flags << "\n"
          << "Contribution Lengths: \n"
          << "GEM Length: " << m_otherContribLen[GEM] << "\n"
          << "OSW Length: " << m_otherContribLen[OSW] << "\n"
          << "AEM Length: " << m_otherContribLen[AEM] << "\n"
          << "Diagnostic Length: " << m_otherContribLen[DIAG] << "\n"
          << "Error Length: " << m_otherContribLen[ERR] << "\n"
          << "TEM Lengths: \n";
        unsigned int i;
        for (i = 0; i < 16; i++) 
            s << "TEM " << i << " Length: " << m_temLen[i] << "\n";
           
        return s;
    }


    inline std::ostream& operator<<(std::ostream &s, const EventSummaryData &obj){
        return obj.fillStream(s);
    }

}//namespace LdfEvent 


#endif
