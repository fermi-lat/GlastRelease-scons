#ifndef RECONEVENT_H
#define RECONEVENT_H

#include "TObject.h"
#include "CalRecon.h"
#include "TkrRecon.h"
#include "AcdRecon.h"

/** @class ReconEvent
 * @brief Primary Root Recon Object - top level object in Root Tree.
 *
 * Portal through which all recon data is accessed.  
 * A ReconEvent object contains:
 * -# Run Number
 * -# Event Id
 * -# A pointer to a TkrRecon object
 * -# A pointer to a CalRecon object
 * -# A pointer to a AcdRecon object 
 *
 * $Header$
*/

class ReconEvent : public TObject
{

public:

    /// Define the bits in the m_eventFlag member
    typedef enum {
        GOOD = 0,
        EVTSEQ = 1
    } EventFlags;

    ReconEvent();

    virtual ~ReconEvent();

    void initialize(UInt_t eventId, UInt_t runId, TkrRecon *tkr, CalRecon *cal, AcdRecon *acd);

    void initEventFlags(UInt_t flags) { m_eventFlags = flags; };

    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    UInt_t getEventId() { return m_eventId; };

    UInt_t getRunId() { return m_runId; };

    AcdRecon* getAcdRecon() { return m_acd; };

    /// provide access to the CAL recon data
    CalRecon* getCalRecon() { return m_cal; };

    ///  Access the TKR recon data
    TkrRecon* getTkrRecon() { return m_tkr; };

    //! Access the Recon Header Flags
   // ReconHeader* getReconFlags() { return &m_recFlags; };
    //! set the ReconHeader pointer 
    //inline void setReconFlags(ReconHeader *r) { m_recFlags = r; };

    UInt_t getEventFlags() const { return m_eventFlags; };
    Bool_t goodEvent() const { return (m_eventFlags == 0); };
    Bool_t badEvent() const { return (m_eventFlags != 0); };
    Bool_t badEventSeq() const { return (m_eventFlags && EVTSEQ); };

private:
	/// Event Id
    UInt_t m_eventId;
	/// Run number
    UInt_t m_runId;
    /// ACD reconstruction data
    AcdRecon *m_acd;
    /// CAL reconstruction data
    CalRecon *m_cal; 
    /// TKR reconstruction data
    TkrRecon *m_tkr;    
    /// pointer to Recon Header data
    //ReconHeader m_recFlags;
    /// Do the easy thing and add flags here
    UInt_t m_eventFlags;

    ClassDef(ReconEvent,3) 
};

#endif

