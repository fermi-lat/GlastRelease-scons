#ifndef SubSystemHeader_H
#define SubSystemHeader_H

#include "TObject.h"

/*! \class SubSystemHeader
\brief ROOT base class to contain the header data for
one event. 
*/
class SubSystemHeader: public TObject {

public:
    enum {
        DeadTimeMask = 0x00002fff,
        DeadTimeCauseMask = 0x001fc000,
        DeadTimeCauseShift = 14
    };

    
    SubSystemHeader();

    virtual ~SubSystemHeader();

    void Clean(Option_t *option="");

    UInt_t getEventId() { return m_eventId; };
    void setEventId(unsigned int id) { m_eventId = id; };

    UInt_t getTimerWord() { return m_timerWord; };
    void setTimerWord(unsigned int time) { m_timerWord = time; };

    UInt_t getTreqStatus() { return m_TREQ_VETO_status; };
    void setTreqStatus(unsigned int stat) { m_TREQ_VETO_status = stat; };

    UInt_t getDeadTime() { return (m_deadtime_and_cause & DeadTimeMask); };

    UChar_t getDeadTimeCause() { return ((m_deadtime_and_cause & DeadTimeCauseMask) >> DeadTimeCauseShift); };

    UInt_t getDeadTimeAndCause() { return m_deadtime_and_cause; };
    void setDeadTimeAndCause(unsigned int t) { m_deadtime_and_cause = t; };

protected:


    UInt_t m_eventId;
    UInt_t m_timerWord;
    UInt_t m_TREQ_VETO_status;
    UInt_t m_deadtime_and_cause;

    ClassDef(SubSystemHeader,1)
};

#endif

