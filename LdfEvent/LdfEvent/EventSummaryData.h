#ifndef LDFCONVERTER_EventSummaryDATA_H
#define LDFCONVERTER_EventSummaryDATA_H

#include <iostream>
#include <vector>

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/IInterface.h"

#include "Event/TopLevel/Definitions.h"
#include "Event/TopLevel/EventModel.h"

#include "eventSummary.h" // from EBF library

#include "enums/EventFlags.h"

/**
* @class EventSummaryData
* @brief TDS for storing the event summary data 
*/

static const CLID& CLID_LdfEventSummaryData = InterfaceID("LdfEventSummaryData", 1, 0);


namespace LdfEvent{

    class EventSummaryData : public DataObject {
    public:

        typedef enum {
            GEM = 0,
            OSW = 1,
            AEM = 2 
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
            for (i = 0; i<16; i++) {
                m_temLen[i] = 0;
                m_diagLen[i] =0;
                m_errLen[i] = 0;
            }
            m_otherContribLen[GEM] = 0;
            m_otherContribLen[AEM] = 0;
            m_otherContribLen[OSW] = 0;
            m_evtSizeInBytes = 0;
        }

        void initialize(unsigned int summary){m_summary=summary;};
        //void initOswEvtSequence(unsigned int evtSeq) { m_evtSequence = evtSeq;};
        void initEventFlags(unsigned int flags) { m_flags = flags; };
        void initTemContribLen(unsigned int *len) {
            unsigned int i;
            for (i = 0; i < 16; i++) { m_temLen[i] = len[i]; }
        }
        void initDiagContribLen(unsigned int *len) {
            unsigned int i;
            for (i = 0; i < 16; i++) { m_diagLen[i] = len[i]; }
        }
        void initErrContribLen(unsigned int *len) {
            unsigned int i;
            for (i = 0; i < 16; i++) { m_errLen[i] = len[i]; }
        }
        void initContribLen(unsigned int *tem, unsigned int gemLen, unsigned int oswLen,
            unsigned int *errLen, unsigned int *diagLen, unsigned int aemLen) {

            m_otherContribLen[GEM] = gemLen;
            m_otherContribLen[OSW] = oswLen;
            m_otherContribLen[AEM] = aemLen;
            initTemContribLen(tem);
            initDiagContribLen(diagLen);
            initErrContribLen(errLen);
        }
  
        void initEventSizeInBytes(unsigned long size) { m_evtSizeInBytes = size; };
        
        unsigned int temLength(unsigned int tem) const { return m_temLen[tem]; }
        unsigned int diagnosticLength(unsigned int tem) const { 
            return m_diagLen[tem]; }
        unsigned int errorLength(unsigned int tem) const { 
            return m_errLen[tem]; }
        unsigned int gemLength() const { return m_otherContribLen[GEM]; }
        unsigned int oswLength() const { return m_otherContribLen[OSW]; }
        unsigned int aemLength() const { return m_otherContribLen[AEM]; }

        const unsigned int* temLength() const { return m_temLen; };
        const unsigned int* diagnosticLength() const { return m_diagLen; };
        const unsigned int* errorLength() const { return m_errLen; };

        void setEvtSeqBit()  { m_flags |= enums::EVTSEQ; };
        void setTkrReconBit() { m_flags |= enums::TKRRECON; };
        unsigned int eventFlags() const { return m_flags; };
        bool goodEvent() const { return (m_flags == 0); };
        bool badEventSeq() const { return ((m_flags & enums::EVTSEQ) != 0); };
        bool badTkrRecon() const { return ((m_flags & enums::TKRRECON) != 0); };
        bool badEvent() const { return (m_flags != 0); } ;
        bool packetError() const { 
            return ((m_flags & enums::PACKETERROR) != 0); };
        bool temError() const { 
            return ((m_flags & enums::SUMMARYERROR) != 0); };
        bool trgParityErrorOr() const { 
            return ((m_flags & enums::TRGPARITYERROR) != 0); };
        bool badLdfStatus() const { 
            return ( (m_flags & enums::BADLDFSTATUS) != 0); };

        //unsigned int eventSequence() const { return m_evtSequence; }
 
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

        unsigned long eventSizeInBytes() const { return m_evtSizeInBytes; };

    private:
    unsigned int m_summary;
    /// Event Flags - used to denote bad events
    /// flags = 0 is "good"
    unsigned int m_flags;
    unsigned int m_temLen[16];
    unsigned int m_diagLen[16];
    unsigned int m_errLen[16];
    unsigned int m_otherContribLen[3];
    // New 32 bit event sequence number from OSW
    //unsigned int m_evtSequence;
    unsigned long m_evtSizeInBytes;
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
          << "TEM Lengths: \n";
        unsigned int i;
        for (i = 0; i < 16; i++) {
            s << "TEM " << i << " Length: " << m_temLen[i] << "\n";
            s << "Diagnostic " << i << " Length: " << m_diagLen[i] << "\n";
            s << "Error " << i << " Length: " << m_errLen[i] << "\n";
        }
           
        return s;
    }


    inline std::ostream& operator<<(std::ostream &s, const EventSummaryData &obj){
        return obj.fillStream(s);
    }

}//namespace LdfEvent 


#endif
