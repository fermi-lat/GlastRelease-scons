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
* - Collection of McParticles
* - Collection of McPosition hits
* - Collection of McIntegrating hits.  
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
    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;
    
    void initialize(UInt_t nEvent, UInt_t nRun);
    
    inline UInt_t getEventId() const { return m_eventId; };
    inline UInt_t getRunId() const { return m_runId; };

    
    /// add a McParticle to list
    void addMcParticle(McParticle *pPart);
    /// retrieve number of McParticles stored in the collection
    inline UInt_t getMcParticleCount() const { 
        return ((m_particleCol) ? m_particleCol->GetEntries() : 0);
    };
    /// get particle from list at index
    McParticle* getMcParticle(UInt_t index) const;
    /// return the full TObjArray containing McParticles
    TObjArray* getMcParticleCol() const { return m_particleCol; };
    
    /// store a new McPositionHit in the collection
    void addMcPositionHit(McPositionHit *hit);
    /// return a McPositionHit corresponding to index
    McPositionHit* getMcPositionHit(UInt_t index) const;
    /// return the number of McPositionHits stored in the collection
    inline UInt_t getMcPositionHitCount() const { 
        return ((m_positionHitCol) ? m_positionHitCol->GetEntries() : 0); 
    };
    /// return the full TObjArray containing McPositionHits
    const TObjArray* getMcPositionHitCol() const { return m_positionHitCol; };

    /// add a new McIntegratingHit to the collection
    void addMcIntegratingHit(McIntegratingHit *hit);
    /// return a McIntegrating hit corresponding to the index 
    McIntegratingHit* getMcIntegratingHit(UInt_t index) const;
    /// return the number of McIntegratingHits stored in the collection
    inline UInt_t getIntegratingHitCount() const { 
        return ((m_integratingHitCol) ? m_integratingHitCol->GetEntries() : 0); 
    };
    /// return the full TObjArray containing McIntegratingHits
    const TObjArray* getMcIntegratingHitCol() const { return m_integratingHitCol; };
    
private:
    /// unique event id for this run
    UInt_t m_eventId;
    
    /// Run number
    UInt_t m_runId;
    
    /// Collection of McParticles
    /// Must be TObjArray due to TRefArray which grows dynamically
    TObjArray *m_particleCol;
    /// static array to allow one-time array creation
    static TObjArray *m_staticParticleCol;
    
    /// Collection of McIntegratingHits
    /// must be a TObjArray due to map in McIntegratingHit
    TObjArray *m_integratingHitCol;
    /// static array to allow one-time array creation
    static TObjArray *m_staticIntegratingHitCol;
    
    /// Collection of McPositionHits
    TObjArray *m_positionHitCol;
    /// static array to allow one-time array creation
    static TObjArray *m_staticPositionHitCol;
    
    ClassDef(McEvent,2) // Monte Carlo Event Class
};

#endif