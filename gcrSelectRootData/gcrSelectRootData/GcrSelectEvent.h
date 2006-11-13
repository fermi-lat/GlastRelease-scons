#ifndef GCRSELECTEVENT_H
#define GCRSELECTEVENT_H

#include "TObject.h"
//#include "CalRecon.h"
//#include "TkrRecon.h"
//#include "AcdRecon.h"
#include "GcrSelect.h"

/** @class GcrSelectEvent
 * @brief Primary Root Recon Object - top level object in Root Tree.
 *
 * Portal through which all recon data is accessed.  
 * A ReconEvent object contains:
 * -# Run Number
 * -# Event Id
 * -# A pointer to a GcrSelect object
 * 	
 * $Header$
*/

class GcrSelectEvent : public TObject
{

public:

    /// Define the bits in the m_eventFlag member
    typedef enum {
        GOOD = 0,
        EVTSEQ = 1
    } EventFlags;

    GcrSelectEvent();

    virtual ~GcrSelectEvent();

    void initialize(UInt_t eventId, UInt_t runId, GcrSelect *gcrSelect);

    void initEventFlags(UInt_t flags) { m_eventFlags = flags; };

    void Clear(Option_t *option="");
    //void Fake( Int_t ievent, Float_t randNum ) ; // for tests
    //Bool_t CompareInRange( const ReconEvent &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    UInt_t getEventId() { return m_eventId; };

    UInt_t getRunId() { return m_runId; };

    //AcdRecon* getAcdRecon() { return m_acd; };

    /// provide access to the CAL recon data
    //CalRecon* getCalRecon() { return m_cal; };

    ///  Access the TKR recon data
    //TkrRecon* getTkrRecon() { return m_tkr; };

    GcrSelect* getGcrSelect() {return m_gcrSelect;}
    

    UInt_t getEventFlags() const { return m_eventFlags; };
    Bool_t goodEvent() const { return (m_eventFlags == 0); };
    Bool_t badEvent() const { return (m_eventFlags != 0); };
    Bool_t badEventSeq() const { return (m_eventFlags && EVTSEQ); };

private:
	/// Event Id
    UInt_t m_eventId;
	/// Run number
    UInt_t m_runId;
    
    GcrSelect *m_gcrSelect;
    /// pointer to Recon Header data
    //ReconHeader m_recFlags;
    /// Do the easy thing and add flags here
    UInt_t m_eventFlags;

    ClassDef(GcrSelectEvent,3) 
};

#endif

