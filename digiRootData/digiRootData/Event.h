
#ifndef Event_H
#define Event_H

//#include "TMath.h"
#include "TObject.h"
#include "L1T.h"
//#include "Tagger.h"

#include "AcdDigi.h"
#include "CalDigi.h"
#include "TkrDigi.h"

/*! \class Event
\brief
 Jan 1999 Daniel Flath - ROOT HTML comments added
 Dec 1999 Daniel Flath - Rewritten for GLAST
 Oct 25,1999 Richard Dubois - Clone from LCD version
 */


class Event: public TObject {
private:
    UInt_t m_run;           // Run number
    UInt_t m_eventId;         // Event number
    
    static CalDigi *m_staticCal;
    static AcdDigi *m_staticAcd;
    static TkrDigi *m_staticTkr;
    static L1T *m_staticL1T;

    CalDigi* m_CAL;       // List of Calorimeter logs
    AcdDigi* m_ACD;       //-> List of ACD tiles
    TkrDigi* m_TKR;       // List of Tracker layers
    static L1T *s_L1T;
    L1T *m_L1Trigger;       // Level 1 trigger object
    //Tagger *m_Tagger;       // Tagger object

public:
    Event();
    virtual ~Event();
    void Clean(Option_t *option="");
    inline void setEventId(UInt_t id) { m_eventId = id; };
    inline void setRun(UInt_t run) { m_run = run; };
    inline void setL1Trigger(L1T *L1TVal) { m_L1Trigger = L1TVal; };
   // inline void setTagger(Tagger *TaggerVal) { m_Tagger = TaggerVal; };
    inline Int_t getEvent() { return m_eventId; };
    inline Int_t getRun() { return m_run; };
    void setCAL(CalDigi *cal) { m_CAL = cal; };
    CalDigi* getCAL() { return m_CAL; };
    void setACD(AcdDigi *acd) { m_ACD = acd; };
    AcdDigi* getACD() { return m_ACD; };
    void setTKR(TkrDigi *tkr){ m_TKR = tkr; };
    TkrDigi* getTKR() { return m_TKR; };
    inline L1T *getL1T() const { return m_L1Trigger; };
    //inline Tagger *getTagger() const { return m_Tagger; };
    

    ClassDef(Event,3)       // Storage for per-event and subsystem data
}; 

#endif





