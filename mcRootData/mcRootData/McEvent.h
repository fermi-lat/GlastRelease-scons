#ifndef MCEVENT_H
#define MCEVENT_H

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "McParticle.h"
#include "McPositionHit.h"
#include "McIntegratingHit.h"


/** @class McEvent
* @brief GLAST Monte Carlo Event class.
*
* The McEvent class contains all the information about a single MC simulated 
* event.  This class contains:
* - Run id
* - Event id
* - List of McParticles
* - List of McPosition hits
* - List of McIntegrating hits.  
* The MC data stored to ROOT files mirrors the data available within the Gaudi 
* TDS.
*  
*  $Header$
*/

class McEvent : public TObject  
{
public:
    
    McEvent();
    
    virtual ~McEvent();
    
    /// clear lists, free pointers, etc., after read from / write to file
    virtual void Clear(Option_t *option ="");
    
    void initialize(Int_t nEvent, Int_t nRun);
    
    inline Int_t getEventId() const { return m_eventId; };
    inline Int_t getRunId() const { return m_eventId; };
    inline Int_t getParticleCount() const { return ((m_particleCol) ? m_particleCol->GetEntries() : 0); };
    
    /// add a McParticle to list
    void addMcParticle(McParticle *pPart);
    /// get particle from list at index
    McParticle* getMcParticle(Int_t index) const;
    
    void addMcPositionHit(McPositionHit *hit);
    McPositionHit* getMcPositionHit(Int_t index) const;

    void addMcIntegratingHit(McIntegratingHit *hit);
    McIntegratingHit* getMcIntegratingHit(Int_t index) const;
    
private:
    /// unique event id for this run
    Int_t m_eventId;
    
    /// Run number
    Int_t m_runId;
    
    /// List of McParticles
    /// Must be TObjArray due to TRefArray which grows dynamically
    TObjArray *m_particleCol;
    /// static array to allow one-time array creation
    static TObjArray *m_staticParticleCol;
    
    /// List of McIntegratingHits
    /// must be a TObjArray due to map in McIntegratingHit
    TObjArray *m_integratingHitCol;
    /// static array to allow one-time array creation
    static TObjArray *m_staticIntegratingHitCol;
    
    /// List of McPositionHits
    TObjArray *m_positionHitCol;
    /// static array to allow one-time array creation
    static TObjArray *m_staticPositionHitCol;
    
    ClassDef(McEvent,2)
};

#endif