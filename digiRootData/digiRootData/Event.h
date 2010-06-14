// Event.h
// Jan 1999 Daniel Flath - ROOT HTML comments added
// Dec 1999 Daniel Flath - Rewritten for GLAST
// Oct 25,1999 Richard Dubois - Clone from LCD version

#ifndef Event_H
#define Event_H

#include "TMath.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TMap.h"
#include "L1T.h"
#include "Tagger.h"
//#include "HouseKeeping.h"

class Event: public TObject {
private:
    UInt_t m_run;           // Run number
    UInt_t m_event;         // Event number
    TObjArray* m_CAL;       // List of Calorimeter logs
    TObjArray* m_ACD;       // List of ACD tiles
    TObjArray* m_TKR;       // List of Tracker layers
    TObjArray* m_MCPart;    // List of MC Particle objects  
    L1T *m_L1Trigger;       // Level 1 trigger object
    Tagger *m_Tagger;       // Tagger object
    //    HouseKeeping *m_HSK;    // HouseKeeping object

public:
    Event();
    ~Event();
    void Clean();
    void Create();
    inline void setEvent(UInt_t event) { m_event = event; };
    inline void setRun(UInt_t run) { m_run = run; };
    inline void setL1Trigger(L1T *L1TVal) { m_L1Trigger = L1TVal; };
    //    inline void setHSK(HouseKeeping *hsk) { m_HSK = hsk; };
    inline void setTagger(Tagger *TaggerVal) { m_Tagger = TaggerVal; };
    inline Int_t getEvent() { return m_event; };
    inline Int_t getRun() { return m_run; };
    TObjArray* getCAL();
    TObjArray* getACD();
    TObjArray* getMCPart();
    TObjArray* getTKR();
    inline L1T *getL1T() const { return m_L1Trigger; };
    inline Tagger *getTagger() const { return m_Tagger; };
    //    inline HouseKeeping *getHSK() const { return m_HSK; };
    //    inline Bool_t getHSKPresent() const { return m_HSK ? kTRUE : kFALSE; };
    

    ClassDef(Event,2)       // Storage for per-event and subsystem data
}; 

#endif





