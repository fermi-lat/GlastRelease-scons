#ifndef SubSystemHeader_H
#define SubSystemHeader_H

#include "TObject.h"

/*! \class SubSystemHeader
\brief Root class to contain the TKR header data for
one event. 
*/
class SubSystemHeader: public TObject {

public:
    
    SubSystemHeader();

    virtual ~SubSystemHeader();

    void Clean(Option_t *option="");

    unsigned long getEventId() { return m_eventId; };
    void setEventId(unsigned int id) { m_eventId = id; };

    unsigned long getTimerWord() { return m_timerWord; };
    void setTimerWord(unsigned int time) { m_timerWord = time; };

    unsigned long getStatus() { return m_TREQ_VETO_status; };
    void setStatus(unsigned int stat) { m_TREQ_VETO_status = stat; };

    unsigned long getDeadTime() { return m_deadtime_and_cause; };

    unsigned long getDeadTimeCause() { return m_deadtime_and_cause; };

    unsigned long getDeadTimeAndCause() { return m_deadtime_and_cause; };
    void setDeadTimeAndCause(unsigned int t) { m_deadtime_and_cause = t; };

protected:

    UInt_t m_eventId;
    UInt_t m_timerWord;
    UInt_t m_TREQ_VETO_status;
    UInt_t m_deadtime_and_cause;

    ClassDef(SubSystemHeader,1)
};

#endif

