#ifndef CalHeader_H
#define CalHeader_H

#include "TObject.h"

/*! \class CalHeader
\brief Root class to contain the CAL header subpacket data for
one event. 
*/
class CalHeader: public TObject {

public:
    
    CalHeader();

    virtual ~CalHeader();

    void Clean(Option_t *option="");

    void setHeader(UInt_t id, UInt_t time, UInt_t stat, UInt_t deadTime, UInt_t cause);

    unsigned long eventId() { return m_eventId; };
    void eventId(unsigned int id) { m_eventId = id; };

    unsigned long timerWord() { return m_timerWord; };
    void timerWord(unsigned int time) { m_timerWord = time; };

    unsigned long status() { return m_TREQ_VETO_status; };
    void status(unsigned int stat) { m_TREQ_VETO_status = stat; };

    unsigned long deadTime() { return m_deadTime; };
    void deadTime(unsigned int time) { m_deadTime = time; };

    unsigned long deadTimeCause() { return m_deadTimeCause; };
    void deadTimeCause(unsigned int cause) { m_deadTimeCause = cause; };

private:

    UInt_t m_eventId;
    UInt_t m_timerWord;
    UInt_t m_TREQ_VETO_status;
    UInt_t m_deadTime;
    UInt_t m_deadTimeCause;
    
    ClassDef(CalHeader,1)
};

#endif

