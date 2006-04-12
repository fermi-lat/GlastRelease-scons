#ifndef ROOT_DigiEvent_H
#define ROOT_DigiEvent_H

#include "TObject.h"
#include "L1T.h"
#include "TClonesArray.h"
#include "TObjArray.h"

#include "AcdDigi.h"
#include "CalDigi.h"
#include "TkrDigi.h"

#include "TkrDiagnosticData.h"
#include "CalDiagnosticData.h"
#include "EventSummaryData.h"
#include "Gem.h"
#include "Tem.h"

#include "MetaEvent.h"
#include "Ccsds.h"

/** @class DigiEvent
 * @brief This is the top-level event class to store the Digi data.
 *
 * The DigiEvent branch contains:
 * - Run Id
 * - Event Id
 * - Time stamp (seconds)
 * - Live Time
 * - Flag denoting if this data was simulated or not
 * - Level One Trigger (L1T)
 * - Collection of AcdDigi objects
 * - Collection of CalDigi objects
 * - Collection of TkrDigi objects
 * - If this is real data, the following are possibly present as well
 *   -# PPC Time Base
 *   -# Original 30Hz clock seconds and nanoseconds
 *   -# EventSummaryData
 *   -# Gem (GLT Electronics Module)
 *   -# CAL Diagnostic trigger primitives from TEMs
 *   -# TKR Diagnostic trigger primitives from TEMs
 *   -# ErrorData accessed throught the TEM collection
 * - The Meta event (aka event context) from FSW, this contains
 *   -# Run information
 *   -# Datagram information
 *   -# Extended context records (aka scalars) with software counters data
 *   -# Event timing information
 *   -# Configuration keys
 * - CCSDS data
 *   -# SCID
 *   -# APID
 *   -# UTC
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
        Double_t liveTime, const L1T& level1, const EventSummaryData &summary, 
        Bool_t fromMc);

    void initialize(UInt_t eventId, UInt_t runId, Double_t time, 
        Double_t liveTime, const L1T& level1, Bool_t fromMc);

    void initEventSummary(const EventSummaryData& summary) { 
         m_summary = summary; };

    void Clear(Option_t *option="");
 
    void Print(Option_t *option="") const;
	
    /** @defgroup DigiEventGroup DigiEvent End-User Interface */
    /*@{*/
    /// Access the DigiEvent event ID 
    inline UInt_t getEventId() const { return m_eventId; };
    /// Access the run number
    inline UInt_t getRunId() const { return m_runId; };

    /// Valid for real or simulated data, in seconds since mission start 
    inline Double_t getTimeStamp() const { return m_timeStamp; };

    /// Valid for both real or simulated data
    inline Double_t getLiveTime() const { return m_liveTime; };

    /// Seconds as reported by original 30 Hz clock (real data only)
    inline UInt_t getEbfTimeSec() const { return m_ebfTimeSec; };
    /// Nanoseconds as reported by original 30 Hz clock (real data only)
    inline UInt_t getEbfTimeNanoSec() const { return m_ebfTimeNanoSec; };

    /// One of two words forming PPC Time Base - see getPpcTimeSeconds
    inline UInt_t getEbfUpperPpcTimeBase() const { return m_ebfUpperPpcTimeBase; };
    /// One of two words forming PPC Time Base - see getPpcTimeSeconds
    inline UInt_t getEbfLowerPpcTimeBase() const { return m_ebfLowerPpcTimeBase; };
    /// Return the approximate number of seconds elapsed since power on
    /// by dividing the value in the PPC registers by 16 MHz
    inline Double_t getEbfPpcTimeSeconds() const { 
        const Double_t sixteenMHz = 16000000.;
        // To handle th 64-bit value - we separate the computation
        // The upper 32 bits would have to be shifted by 31 (or multiplied)
        // by 2^32 - we divide this by 16000000 to get the upperMultiplier
        const Double_t upperMultiplier = 268.435456; 
        Double_t lower = m_ebfLowerPpcTimeBase / sixteenMHz;
        Double_t upper = m_ebfUpperPpcTimeBase * upperMultiplier;
        return (upper + lower);
    };

    /// New name for getEbfPpcTimeSeconds()
    inline Double_t getPpcTimeSeconds() const { return getEbfPpcTimeSeconds(); };

    /// Flag denoting if this event was generated from a Monte Carlo Simulation
    inline Bool_t getFromMc() const { return m_fromMc; };

    /// retrieve the whole TClonesArray of Acd Digi data
    const TClonesArray* getAcdDigiCol() const { return m_acdDigiCol; };
    /// retrieve a specific AcdDigi, based upon the layer, face, row, col, if not found returns null
    const AcdDigi* getAcdDigi(UShort_t l, UShort_t f, UShort_t r, UShort_t c) const;
    /// retrieve a specific AcdTile, based upon a valid AcdId, if not found returns null
    const AcdDigi* getAcdDigi(const AcdId &id) const;

    /// retrieve the whole CalDigi collection
    const TClonesArray* getCalDigiCol() const;
    /// retrieve one CalDigi object from the collection, using the index into the array
    const CalDigi* getCalDigi(UInt_t i) const;

    /// retrieve the whole TObjArray of TkrDigi Data
    const TObjArray* getTkrDigiCol() const { return m_tkrDigiCol; };
    /// retrieve a TkrDigi from the collection, using the index into the array
    const TkrDigi* getTkrDigi(UInt_t i) const;

    /// Access Level 1 Trigger data const
    inline const L1T& getL1T() const { return m_levelOneTrigger; };    
    /// Access Level 1 Trigger data non-const
    inline L1T& getL1T() { return m_levelOneTrigger; };    

    /// Returns a reference to the EventSummaryData (valid real data only) const
    inline const EventSummaryData& getEventSummaryData() const { return m_summary; };
    /// Returns a reference to the EventSummaryData (valid real data only) non-const
    inline EventSummaryData& getEventSummaryData() { return m_summary; };

    /// Returns the whole CalDiagnostic Collection
    const TClonesArray *getCalDiagnosticCol() const { return m_calDiagnosticCloneCol;};
    /// Returns one CalDiagnostic object located at index i in the collection
    const CalDiagnosticData* getCalDiagnostic(UInt_t i) const;

    /// Returns the whole TkrDiagnostic Collection
    const TClonesArray *getTkrDiagnosticCol() const { return m_tkrDiagnosticCloneCol;};
    /// Returns one TkrDiagnostic object located at index i in the collection
    const TkrDiagnosticData* getTkrDiagnostic(UInt_t i) const;

    /// Returns a reference to the Gem
    const Gem& getGem() const { return m_gem; };

    /// Return one Tem object located at index i in the collection
    const Tem* getTem(Int_t i) const;

    /// Returns whole Tem Collection
    const TClonesArray *getTemCol() const { return m_temCloneCol;};
	/*@}*/

    /// Returns the MetaEvent
    const MetaEvent& getMetaEvent() const { return m_metaEvent; }

    /// Copies in the MetaEvent
    void setMetaEvent(const MetaEvent& meta) { m_metaEvent = meta; }

    const Ccsds& getCcsds() const { return m_ccsds; }

    void setCcsds(const Ccsds& c) { m_ccsds = c; }

    /// clear the whole array (necessary because of the consts-s)
    void clearTkrDigiCol() { m_tkrDigiCol->Clear(); }
    /// Add a TkrDigi into the TKR data collection
    void addTkrDigi(TkrDigi *digi);

    CalDigi* addCalDigi();
	
    /// Add a new AcdDigi entry into the ACD digi array
    AcdDigi* addAcdDigi(const AcdId& id, const VolumeIdentifier& volId, 
        Float_t energy, UShort_t *pha, Bool_t *veto, Bool_t *low, Bool_t *high);

    void initGem(const Gem& gem) { m_gem = gem; };
    TkrDiagnosticData* addTkrDiagnostic();
    CalDiagnosticData* addCalDiagnostic();

    Tem* addTem();

    void setEbfTime(UInt_t timeSec, UInt_t timeNanoSec, 
                    UInt_t upperPpcTimeBase, UInt_t lowerPpcTimeBase) {
        m_ebfTimeSec = timeSec;
        m_ebfTimeNanoSec = timeNanoSec;
        m_ebfUpperPpcTimeBase = upperPpcTimeBase;
        m_ebfLowerPpcTimeBase = lowerPpcTimeBase;
    }


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
    // Avoid re-use of m_calDigiCol name that used to be a TObjArray*
    TObjArray *m_calDigiCol; 
    TClonesArray *m_calDigiCloneCol; //->
    Int_t m_numCalDigis;
    static TClonesArray *s_calDigiStaticCol; //! Collection of Cal Digi objects

    /// data members to store TKR data
    static TObjArray *s_staticTkrDigiCol;
    TObjArray* m_tkrDigiCol; //-> List of Tracker layers


    EventSummaryData m_summary;

    TClonesArray *m_tkrDiagnosticCloneCol; //->
    Int_t m_numTkrDiagnostics;
    static TClonesArray *s_tkrDiagnosticStaticCol; //! Collection of TKR diagnostic for EM

    TClonesArray *m_calDiagnosticCloneCol; //->
    Int_t m_numCalDiagnostics;
    static TClonesArray *s_calDiagnosticStaticCol; //! Collection of CAL dianostics for EM

    /// These two words represent the original clock for the EM at 30 Hz
    UInt_t m_ebfTimeSec;
    UInt_t m_ebfTimeNanoSec;


    /// These two words represent the contents of 2 32-bit PPC registers
    /// Concatenate these two words into one 64 bit word divide by 16 MHz
    /// and get an approximate time to be used to determine the sequence of 
    /// events
    UInt_t m_ebfUpperPpcTimeBase;
    UInt_t m_ebfLowerPpcTimeBase;

    Gem m_gem;
    Double_t m_liveTime;

    TClonesArray *m_temCloneCol; //->
    Int_t m_numTem;
    static TClonesArray *s_temStaticCol; //!

    MetaEvent m_metaEvent; 

    Ccsds m_ccsds;

    ClassDef(DigiEvent,16) // Storage for Raw(Digi) event and subsystem data
}; 
 
#endif





