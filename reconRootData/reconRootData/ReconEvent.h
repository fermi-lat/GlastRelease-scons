#ifndef RECON_H
#define RECON_H

#include "TObject.h"
#include "CalRecon.h"
#include "TkrRecon.h"
#include "AcdRecon.h"

/** @class Recon
 * @brief Primary Root Recon Object - top level object in Root Tree
 * Top level Root Recon object - portal through which all recon data is
   accessed.  
   A Recon object contains 3 pointers:
   1) A pointer to a TkrRecon object
   2) A pointer to a CalRecon object
   3) A pointer to a ReconHeader object 
*/

class ReconEvent : public TObject
{

public:

    ReconEvent();

    virtual ~ReconEvent();

    void initialize(UInt_t eventId, UInt_t runId, TkrRecon *tkr, CalRecon *cal, AcdRecon *acd);

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

private:
    UInt_t m_eventId;
    UInt_t m_runId;
    /// ACD reconstruction data
    AcdRecon *m_acd;
    /// CAL reconstruction data
    CalRecon *m_cal; 
    /// TKR reconstruction data
    TkrRecon *m_tkr;    
    /// pointer to Recon Header data
    //ReconHeader m_recFlags;

    ClassDef(ReconEvent,2) 
};

#endif

