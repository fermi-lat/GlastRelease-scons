// RMcVertex.h
//
// The RMcEvent class contains all the information about a single MC simulated event
// This includes run/event numbers, and a list of McParticles generated during the event.
// Developed initially for the GLAST Ballon Flight sim.
//
// May 30 2001 Daniel Flath - Creation

#ifndef RMCEVENT_H
#define RMCEVENT_H

#include "TObject.h"
#include "TObjArray.h"

#include "RMcVertex.h"


/*! GLAST Monte Carlo vertex class.  Designed originally for 2001 balloon flight.
 *  Gives position information for McParticle which keeps a pointer to a vertex.
 *  
 *  \author Daniel Flath
 *  
 *  \b Revisions:
 *  - 30 May 2001   Daniel Flath    Creation
 *  - 06 Jun 2001   Daniel Flath    Final rewrite for integration with ROOTWriter
 *  - 14 Jun 2001   Daniel Flath    DOxygen style comments updated for checkin with ROOTWriter
 */

class RMcEvent : public TObject  
{
private:
    //! Event number (unique id)
    Int_t m_nEvent;

    //! Run number
    Int_t m_nRun;

    //! List of particles
    /*! Ensure that if list is modified, that m_nPartCount is updated accordingly
     */
    TObjArray *m_pVertList;

    //! Number of particles in list
    /*! Must be updated to reflect modifications to the m_pPartList
     */
    Int_t m_nVertCount;
public:
    ////////////////////// construction/destruction: ///////////////////////

    //! default constructor. OK TO USE -- see details.
    RMcEvent();

    //! destructor
    ~RMcEvent();

    //////////////////////   data access functions:  ///////////////////////

    //! clear lists, free pointers, etc., after read from / write to file
    void Clean();

    //! recreate lists, nullify pointers, etc., before read from / write to file
    void Create();

    //! get event number
    inline Int_t getEvent()     const { return m_nEvent; };

    //! get run number
    inline Int_t getRun()       const { return m_nEvent; };

    //! get number of particles in list
    inline Int_t getVertCount() const { return m_nVertCount; };

    //! set event number
    inline void setEvent(Int_t nEvent)  { m_nEvent = nEvent; };

    //! set run number
    inline void setRun(Int_t nRun)      { m_nRun = nRun; };

    //! get particle from list at index
    RMcVertex *getVert(Int_t index) const;

    //! add a particle to list
    Int_t addVert(RMcVertex *pPart);


    ClassDef(RMcEvent,1)
};

#endif