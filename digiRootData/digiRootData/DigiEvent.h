#ifndef ROOT_DigiEvent_H
#define ROOT_DigiEvent_H

#include "TObject.h"
#include "L1T.h"
#include "TClonesArray.h"
#include "TObjArray.h"

#include "AcdDigi.h"
#include "CalDigi.h"
#include "TkrDigi.h"

/** @class DigiEvent
 * @brief This is the top-level event class to store the Digi data.
 *
 * The DigiEvent branch contains:
 * - Run Id
 * - Event Id
 * - Time stamp
 * - Flag denoting if this data was simulated or not
 * - Level One Trigger
 * - Collection of AcdDigi objects
 * - Collection of CalDigi objects
 * - Collection of TkrDigi objects
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

    void initialize(UInt_t eventId, UInt_t runId, Double_t time, 
        const L1T& level1, Bool_t fromMc=true);

    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    inline Double_t getTimeStamp() { return m_timeStamp; };

    /// Access the DigiEvent number
    inline UInt_t getEventId() { return m_eventId; };

    /// Access the run number
    inline UInt_t getRunId() { return m_runId; };

    inline Bool_t getFromMc() { return m_fromMc; };

    /// retrieve the whole TClonesArray of Acd Digi data
    const TClonesArray* getAcdDigiCol() const { return m_acdDigiCol; };
    /// Add a new AcdDigi entry into the ACD digi array
    AcdDigi* addAcdDigi(const AcdId& id, const VolumeIdentifier& volId, 
        Float_t energy, UShort_t *pha, Bool_t *veto, Bool_t *low, Bool_t *high);

    /// retrieve a specific AcdDigi, based upon the layer, face, row, col, if not found returns null
    const AcdDigi* getAcdDigi(UShort_t l, UShort_t f, UShort_t r, UShort_t c) const;

    /// retrieve a specific AcdTile, based upon a valid AcdId, if not found returns null
    const AcdDigi* getAcdDigi(const AcdId &id) const;

    /// retrieve the CalDigi object
    const TClonesArray* getCalDigiCol() { return m_calDigiCol; };
    CalDigi* addCalDigi();
    const CalDigi* getCalDigi(UInt_t i) const;

    /// retrieve the whole TObjArray of TkrDigi Data
    const TObjArray* getTkrDigiCol() const { return m_tkrDigiCol; };
    /// Add a TkrDigi into the TKR data collection
    void addTkrDigi(TkrDigi *digi);
    /// retrieve a TkrDigi from the collection, using the index into the array
    const TkrDigi* getTkrDigi(UInt_t i) const;

    /// Access Level 1 Trigger data
    inline const L1T& getL1T() const { return m_levelOneTrigger; };    

private:
    /// Time in seconds
    Double_t m_timeStamp;

    /// Run number
    UInt_t m_runId;  
    /// Event Number 
    UInt_t m_eventId;  
    
    /// Denote whether or not this data was simulated
    Bool_t m_fromMc;

    /// Store Level 1 trigger
    L1T m_levelOneTrigger;

    /// data members to store ACD data
    TClonesArray *m_acdDigiCol;  //-> 
    Int_t m_numAcdDigis;
    static TClonesArray *s_acdDigiStaticCol; //!

    /// data members to store CAL data
    TClonesArray *m_calDigiCol; //->
    Int_t m_numCalDigis;
    static TClonesArray *s_calDigiStaticCol; //! Collection of Cal Digi objects

    /// data members to store TKR data
    static TObjArray *s_staticTkrDigiCol;
    TObjArray* m_tkrDigiCol; //-> List of Tracker layers

    ClassDef(DigiEvent,6) // Storage for Raw(Digi) event and subsystem data
}; 

#endif





