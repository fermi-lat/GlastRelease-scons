
#ifndef Event_H
#define Event_H

#include "TObject.h"
#include "L1T.h"
//#include "Tagger.h"
#include "TClonesArray.h"
#include "TObjArray.h"

#include "AcdHeader.h"
#include "CalHeader.h"
#include "TkrHeader.h"

#include "AcdTile.h"
#include "CalLog.h"
#include "TkrLayer.h"

/*! \class Event
\brief This is the top-level event class to store the
raw (digi) data.
 Jun 2001 Heather Kelly - revised to use TClonesArray
 Jan 2000 Daniel Flath - ROOT HTML comments added
 Dec 1999 Daniel Flath - Rewritten for GLAST
 Oct 25,1999 Richard Dubois - Clone from LCD version
 */


class Event: public TObject {
private:
    /// Run number
    UInt_t m_run;  
    /// Event Number 
    UInt_t m_eventId;  
    
    /// Store Level 1 trigger
    L1T m_L1T;

    //Tagger *m_Tagger;       // Tagger object

    /// data members to store ACD data
    TClonesArray *m_AcdDigiVec;  //-> 
    Int_t m_numTiles;
    static TClonesArray *m_staticAcdDigiVec;
    AcdHeader m_AcdHeader;

    /// data members to store XGT data
    TClonesArray *m_XgtDigiVec;  //-> 
    Int_t m_numXgts;
    static TClonesArray *m_staticXgtDigiVec;

    /// data members to store CAL data
    TClonesArray *m_CalDigiVec;  //->
    Int_t m_numLogs;
    static TClonesArray *m_staticCalDigiVec;
    CalHeader m_CalHeader;

    /// data members to store TKR data
    static TObjArray *m_staticTkrDigiVec;
    TObjArray* m_TkrDigiVec;       //-> List of Tracker layers
    Int_t m_numLayers;
    TkrHeader m_TkrHeader;


public:

    Event();
    virtual ~Event();

    void Clean(Option_t *option="");
    

    /// Access the Event number
    inline void setEventId(UInt_t id) { m_eventId = id; };
    inline Int_t getEventId() { return m_eventId; };

    /// Access the run number
    inline void setRunId(UInt_t run) { m_run = run; };
    inline Int_t getRunId() { return m_run; };

    /// Access ACD data
    AcdHeader* getAcdHeader() { return &m_AcdHeader; };
    /// retrieve the whole TClonesArray of Acd Digi data
    TClonesArray* getAcdDigi() { return m_AcdDigiVec; };
    /// Add a new AcdTile entry into the ACD digi array
    AcdTile* addAcdTile(Short_t id);
    /// retrieve a specific AcdTile - if not found, returns null
    const AcdTile* getAcdTile(UShort_t id);

    /// retrieve the whole TClonesArray of XGT Digi data
    TClonesArray* getXgtDigi() { return m_XgtDigiVec; };
    /// add a new AcdTile to store XGT data into the XGT digi array
    AcdTile* addXgt(Short_t id);
    /// retrieve a specific XGT (AcdTile), if not found, returns null
    const AcdTile* getXgt(UShort_t id);

    /// Access CAL data
    CalHeader* getCalHeader() { return &m_CalHeader; };
    /// retrieve the whole TClonesArray of Cal Digi data
    TClonesArray* getCalDigi() { return m_CalDigiVec; };
    /// Add a new CalLog entry into the CAL digi array
    CalLog* addCalLog();

    /// Access TKR data
    TkrHeader* getTkrHeader() { return &m_TkrHeader; };
    /// retrieve the whole TObjArray of Tkr Digi Data
    TObjArray* getTkrDigi() { return m_TkrDigiVec; };
    /// Add a TkrLayer entry into the TKR data array
    void addTkrLayer(TkrLayer *layer);
    /// retrieve a specific layer of data, identified by Layer Number
    const TkrLayer* getTkrLayer(unsigned int layerNum);

    /// Access Level 1 Trigger data
    inline L1T* getL1T() { return &m_L1T; };

    // inline void setTagger(Tagger *TaggerVal) { m_Tagger = TaggerVal; };
    //inline Tagger *getTagger() const { return m_Tagger; };
    
    ClassDef(Event,3)       // Storage for Raw(Digi) event and subsystem data
}; 

#endif





