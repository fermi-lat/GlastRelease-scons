
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
#include "TkrDigi.h"
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


public:

    DigiEvent();
    virtual ~DigiEvent();

    void initialize(UInt_t eventId, UInt_t runId, Bool_t fromMc=true);

    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    /// Access the DigiEvent number
    inline UInt_t getEventId() { return m_eventId; };

    /// Access the run number
    inline UInt_t getRunId() { return m_runId; };

    inline Bool_t getFromMc() { return m_fromMc; };

    /// Access ACD data
    //AcdHeader* getAcdHeader() { return &m_AcdHeader; };
    /// retrieve the whole TClonesArray of Acd Digi data
    //TClonesArray* getAcdDigi() { return m_AcdDigiVec; };
    /// Add a new AcdTile entry into the ACD digi array
    //AcdTile* addAcdTile(UInt_t id, short base=10, short used=1);
    /// retrieve a specific AcdTile - if not found, returns null
    /// User should provide a valid (Ritz) Acd id, in base 10 - where 010 == 10
    //const AcdTile* getAcdTile(UInt_t id);
    /// retrieve a specific AcdTile, based upon the layer, face, row, col, if not found returns null
    //const AcdTile* getAcdTile(short l, short f, short r, short c);
    /// retrieve a specific AcdTile, based upon a valid AcdId, if not found returns null
    //const AcdTile* getAcdTile(AcdId &id);

    /// Access CAL data
    //CalHeader* getCalHeader() { return &m_CalHeader; };
    /// retrieve the CalDigi object
    const TObjArray* getCalDigiCol() { return m_calDigiCol; };
    void addCalDigi(CalDigi *cal);
    const CalDigi* getCalDigi(UInt_t i) const;

    /// Access TKR data
    //TkrHeader* getTkrHeader() { return &m_TkrHeader; };
    /// retrieve the whole TObjArray of Tkr Digi Data
    TObjArray* getTkrDigiCol() { return m_tkrDigiCol; };
    /// Add a TkrDigi into the TKR data collection
    void addTkrDigi(TkrDigi *digi);
    /// retrieve a TkrDigi from the collection, using the index into the array
    const TkrDigi* getTkrDigi(UInt_t i) const;

    /// Access Level 1 Trigger data
    //inline L1T* getL1T() { return &m_L1T; };

    //inline LiveTime* getLiveTime() { return &m_liveTime; };
    

private:
    /// Run number
    UInt_t m_runId;  
    /// Event Number 
    UInt_t m_eventId;  
    
    /// Denote whether or not this data was simulated
    Bool_t m_fromMc;

    /// Store Level 1 trigger
    //L1T m_L1T;

    /// data members to store ACD data
    //TClonesArray *m_AcdDigiVec;  //-> 
    //Int_t m_numTiles;
    //static TClonesArray *m_staticAcdDigiVec;
    //AcdHeader m_AcdHeader;

    /// data members to store CAL data
    TObjArray *m_calDigiCol; //->
    static TObjArray *s_calDigiStaticCol; //! Collection of Cal Digi objects
    //CalHeader m_CalHeader;

    /// data members to store TKR data
    static TObjArray *s_staticTkrDigiCol;
    TObjArray* m_tkrDigiCol; //-> List of Tracker layers
    //Int_t m_numLayers;
    //TkrHeader m_TkrHeader;

    /// store LiveTime counter data
    //LiveTime m_liveTime;
    ClassDef(DigiEvent,4) // Storage for Raw(Digi) event and subsystem data
}; 

#endif





