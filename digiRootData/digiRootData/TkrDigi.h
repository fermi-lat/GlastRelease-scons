#ifndef TkrDigi_H
#define TkrDigi_H

#include "TObject.h"
#include "TObjArray.h"

#include "TkrLayer.h"

/*! \class TkrDigi
\brief Root class to contain the TKR digitization data for
one event. 
*/
class TkrDigi: public TObject {

public:
    
    TkrDigi();

    virtual ~TkrDigi();

    void Clean(Option_t *option="");

    /// provide mechanism to update array of logs
    void Add(TkrLayer *layer);

    /// provide access to array of logs
    TObjArray* getLayers() { return m_layers; };

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

    TObjArray *m_layers;

    static TObjArray *m_staticArray;

    unsigned long m_eventId;
    unsigned long m_timerWord;
    unsigned long m_TREQ_VETO_status;
    unsigned long m_deadTime;
    unsigned long m_deadTimeCause;
    
    ClassDef(TkrDigi,1)
};

#endif

