#ifndef CalDigi_H
#define CalDigi_H

#include "TObject.h"
#include "TClonesArray.h"

#include "CalHit.h"
#include "LogID.h"

/*! \class CalDigi
\brief Root class to contain the CAL digitization data for
one event. 
*/
class CalDigi: public TObject {

public:
    
    CalDigi();

    virtual ~CalDigi();

    void Clean(Option_t *option="");

    /// provide mechanism to update array of logs
    CalHit* Add(LogID *id);

    /// provide access to array of logs
    TClonesArray* getLogs() { return m_logs; };

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

    TClonesArray *m_logs; //->

    static TClonesArray *m_staticCalArray;

    unsigned int m_numLogs;

    unsigned long m_eventId;
    unsigned long m_timerWord;
    unsigned long m_TREQ_VETO_status;
    unsigned long m_deadTime;
    unsigned long m_deadTimeCause;
    
    ClassDef(CalDigi,1)
};

#endif

