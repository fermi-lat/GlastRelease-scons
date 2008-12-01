#ifndef overlayRootData_EventOverlay_H
#define overlayRootData_EventOverlay_H

#include "TObject.h"
#include "TClonesArray.h"
#include "TObjArray.h"

#include "AcdOverlay.h"
#include "CalOverlay.h"
#include "TkrOverlay.h"
#include "GemOverlay.h"

/** @class EventOverlay
 * @brief This is the top-level event class to store the Overlay data.
 *
 * $Header$
 */

class EventOverlay: public TObject 
{
public:

    EventOverlay();
    virtual ~EventOverlay();

    void initialize(UInt_t eventId, UInt_t runId, Double_t time, 
        Double_t liveTime, Bool_t fromMc);

    void setGemOverlay(const GemOverlay& gem)  {m_gemOverlay = gem;}

    void Clear(Option_t *option="");
 
    void Print(Option_t *option="") const;
    
    /** @defgroup EventOverlayGroup EventOverlay End-User Interface */
    /*@{*/
    /// Access the EventOverlay event ID 
    inline UInt_t     getEventId() const       { return m_eventId; };
    /// Access the run number
    inline UInt_t     getRunId() const         { return m_runId; };

    /// Valid for real or simulated data, in seconds since mission start 
    inline Double_t   getTimeStamp() const     { return m_timeStamp; };

    /// Currently valid only for simulated data
    inline Double_t   getLiveTime() const      { return m_liveTime; };

    /// Flag denoting if this event was generated from a Monte Carlo Simulation
    inline Bool_t     getFromMc() const        { return m_fromMc; };

    /// Returns a reference to the Gem
    const GemOverlay& getGemOverlay() const    { return m_gemOverlay; };

    /// retrieve the whole TClonesArray of Acd Digi data
    const TObjArray*  getAcdOverlayCol() const { return m_acdOverlayCol; };
    /// retrieve one AcdDigi object from the collection, using the index into the array
    const AcdOverlay* getAcdOverlay(UInt_t i) const;
    /// retrieve the whole CalDigi collection
    const TObjArray*  getCalOverlayCol() const { return m_calOverlayCol;};
    /// retrieve one CalDigi object from the collection, using the index into the array
    const CalOverlay* getCalOverlay(UInt_t i) const;
    /// retrieve the whole TObjArray of TkrDigi Data
    const TObjArray*  getTkrOverlayCol() const { return m_tkrOverlayCol; };
    /// retrieve a TkrDigi from the collection, using the index into the array
    const TkrOverlay* getTkrOverlay(UInt_t i) const;

    /// clear the whole array (necessary because of the consts-s)
    void clearTkrOverlayCol() { m_tkrOverlayCol->Clear(); }
    /// Add a TkrOverlay into the TKR data collection
    void addTkrOverlay(TkrOverlay *tkrOverlay);

    /// clear the whole array (necessary because of the consts-s)
    void clearCalOverlayCol() { m_calOverlayCol->Clear(); }
    /// Add a CalOverlay into the CAL data collection
    void addCalOverlay(CalOverlay* calOverlay);

    /// clear the whole array (necessary because of the consts-s)
    void clearAcdOverlayCol() { m_acdOverlayCol->Clear(); }
    /// Add a AcdOverlay into the ACD data collection
    void addAcdOverlay(AcdOverlay* acdOverlay);

private:
    /// Time in seconds
    Double_t   m_timeStamp;

    /// Run number
    UInt_t     m_runId;  
    /// Event Number 
    UInt_t     m_eventId;  
    
    /// Denote whether or not this data was simulated
    Bool_t     m_fromMc;

    /// livetime
    Double_t   m_liveTime;

    /// data members to store ACD data
    TObjArray* m_acdOverlayCol;  
    Int_t      m_numAcdOverlays;

    /// data members to store CAL data
    // Avoid re-use of m_calDigiCol name that used to be a TObjArray*
    TObjArray* m_calOverlayCol; 
    Int_t      m_numCalOverlays;

    /// data members to store TKR data
    TObjArray* m_tkrOverlayCol; //-> List of Tracker layers
    Int_t      m_numTkrOverlays;

    /// Storage for gem information
    GemOverlay m_gemOverlay;

    ClassDef(EventOverlay,1) // Storage for Raw(Digi) event and subsystem data
}; 
 
#endif





