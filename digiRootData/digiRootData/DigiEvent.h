
#ifndef DigiEvent_H
#define DigiEvent_H

#include "TObject.h"
#include "L1T.h"
#include "TClonesArray.h"
#include "TObjArray.h"

#include "AcdHeader.h"
#include "CalHeader.h"
#include "TkrHeader.h"

#include "AcdTile.h"
#include "CalDigi.h"
#include "TkrLayer.h"

#include "LiveTime.h"

/** @class DigiEvent
 * @brief This is the top-level event class to store the raw (digi) data.
 *
 * @li Jun 2001 Heather Kelly - revised to use TClonesArray
 * @li Jan 2000 Daniel Flath - ROOT HTML comments added
 * @li Dec 1999 Daniel Flath - Rewritten for GLAST
 * @li Oct 25,1999 Richard Dubois - Clone from LCD version
 * $Header$
 */

class DigiEvent: public TObject {
private:
    /// Run number
    UInt_t m_run;  
    /// Event Number 
    UInt_t m_eventId;  
    
    /// Store Level 1 trigger
    L1T m_L1T;

    /// data members to store ACD data
    TClonesArray *m_AcdDigiVec;  //-> 
    Int_t m_numTiles;
    static TClonesArray *m_staticAcdDigiVec;
    AcdHeader m_AcdHeader;

    /// data members to store CAL data
    CalDigi *m_calDigi;
    CalHeader m_CalHeader;

    /// data members to store TKR data
    static TObjArray *m_staticTkrDigiVec;
    TObjArray* m_TkrDigiVec;       //-> List of Tracker layers
    Int_t m_numLayers;
    TkrHeader m_TkrHeader;

    /// store LiveTime counter data
    LiveTime m_liveTime;

public:

    DigiEvent();
    virtual ~DigiEvent();

    void Clear(Option_t *option="");

    void Print(Option_t *option="");

    /// Access the DigiEvent number
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
    AcdTile* addAcdTile(UInt_t id, short base=10, short used=1);
    /// retrieve a specific AcdTile - if not found, returns null
    /// User should provide a valid (Ritz) Acd id, in base 10 - where 010 == 10
    const AcdTile* getAcdTile(UInt_t id);
    /// retrieve a specific AcdTile, based upon the layer, face, row, col, if not found returns null
    const AcdTile* getAcdTile(short l, short f, short r, short c);
    /// retrieve a specific AcdTile, based upon a valid AcdId, if not found returns null
    const AcdTile* getAcdTile(AcdId &id);

    /// Access CAL data
    CalHeader* getCalHeader() { return &m_CalHeader; };
    /// retrieve the CalDigi object
    const CalDigi* getCalDigi() { return m_calDigi; };

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

    inline LiveTime* getLiveTime() { return &m_liveTime; };
    
    ClassDef(DigiEvent,4) // Storage for Raw(Digi) event and subsystem data
}; 

#endif





