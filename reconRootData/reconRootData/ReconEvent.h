#ifndef RECON_H
#define RECON_H

#include "TObject.h"
#include "CalRecon.h"

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

    void initialize(UInt_t eventId, UInt_t runId);

    void Clear(Option_t *option="");

    void Print(Option_t *option="");

    UInt_t getEventId() { return m_eventId; };

    UInt_t getRunId() { return m_runId; };

    //! provide access to the CAL recon data
    const CalRecon* getCalRecon() const { return &m_cal; };
    CalRecon* getCalRecon() { return &m_cal; };
    //! set the CalRecon pointer
    //inline void setCalRecon(CalRecon *cal) { m_cal = cal; };

    //!  Access the TKR recon data
    //TkrRecon* getTkrRecon() { return &m_tkr; };
    //! Set the TkrRecon pointer
    //inline void setTkrRecon(TkrRecon *tkr) { m_tkr = tkr; };

    //! Access the Recon Header Flags
   // ReconHeader* getReconFlags() { return &m_recFlags; };
    //! set the ReconHeader pointer 
    //inline void setReconFlags(ReconHeader *r) { m_recFlags = r; };

private:
    UInt_t m_eventId;
    UInt_t m_runId;
    /// CAL reconstruction data
    CalRecon m_cal;    
    /// pointer to TKR reconstruction data
   // TkrRecon m_tkr;    
    /// pointer to Recon Header data
    //ReconHeader m_recFlags;

    ClassDef(ReconEvent,1) 
};

#endif

