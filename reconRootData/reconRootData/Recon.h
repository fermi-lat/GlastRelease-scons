#ifndef RECON_H
#define RECON_H

// ROOT includes
#include "TObject.h"
#include "CalRecon.h"
#include "TkrRecon.h"
#include "ReconHeader.h"

//! Primary Root Recon Object - top level object in Root Tree
/*! Top level Root Recon object - portal through which all recon data is
    accessed.  
    A Recon object contains 3 pointers:
    1) A pointer to a TkrRecon object
    2) A pointer to a CalRecon object
    3) A pointer to a ReconHeader object */

class Recon : public TObject
{
private:
    //! pointer to CAL reconstruction data
    CalRecon m_cal;    
    //! pointer to TKR reconstruction data
    TkrRecon m_tkr;    
    //! pointer to Recon Header data
    ReconHeader m_recFlags;

public:

    //! constructor
    Recon();

    //! destructor
    virtual ~Recon();

    //! Deletes all dynamically allocated memory
    void Clean();
    //! Creates the TkrRecon and CalRecon objects
    void Create();

    //! provide access to the CAL recon data
    CalRecon* getCalRecon() { return &m_cal; };
    //! set the CalRecon pointer
    //inline void setCalRecon(CalRecon *cal) { m_cal = cal; };

    //!  Access the TKR recon data
    TkrRecon* getTkrRecon() { return &m_tkr; };
    //! Set the TkrRecon pointer
    //inline void setTkrRecon(TkrRecon *tkr) { m_tkr = tkr; };

    //! Access the Recon Header Flags
    ReconHeader* getReconFlags() { return &m_recFlags; };
    //! set the ReconHeader pointer 
    //inline void setReconFlags(ReconHeader *r) { m_recFlags = r; };


    ClassDef(Recon,2) 
};

#endif

