#ifndef TkrHeader_H
#define TkrHeader_H

#include "TObject.h"

/*! \class TkrHeader
\brief Root class to contain the TKR header data for
one event. 
*/
class TkrHeader: public TObject {

public:
    
    TkrHeader();

    virtual ~TkrHeader();

    void Clean(Option_t *option="");

    UInt_t eventId() { return m_eventId; };
    void eventId(UInt_t id) { m_eventId = id; };

    UInt_t timerWord() { return m_timerWord; };
    void timerWord(UInt_t time) { m_timerWord = time; };

    UInt_t status() { return m_TREQ_VETO_status; };
    void status(UInt_t stat) { m_TREQ_VETO_status = stat; };

    UInt_t deadTime() { return m_deadTime; };
    void deadTime(UInt_t time) { m_deadTime = time; };

    UInt_t deadTimeCause() { return m_deadTimeCause; };
    void deadTimeCause(UInt_t cause) { m_deadTimeCause = cause; };

private:

    UInt_t m_eventId;
    UInt_t m_timerWord;
    UInt_t m_TREQ_VETO_status;
    UInt_t m_deadTime;
    UInt_t m_deadTimeCause;
    
    ClassDef(TkrHeader,1)
};

#endif

