#ifndef MCEVENT_H
#define MCEVENT_H

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "McParticle.h"
#include "McPositionHit.h"
#include "McIntegratingHit.h"
#include "McTrajectory.h"

/** @class McEvent
* @brief GLAST Monte Carlo Event class.
*
* The McEvent class contains all the information about a single MC simulated 
* event.  This class contains:
* - Run id
* - Event id
* - Collection of McParticle
* - Collection of McPositionHit
* - Collection of McIntegratingHit
* - Collection of McTrajectories
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
    // DC: I CANNOT MAKE THE USUAL Fake() and CompareInRange(),
    // because McEvent is designed as a singleton.
    void Fake( Int_t ievent, Int_t irun, Float_t randNum ) ; // for tests
    Bool_t CompareToFake( Int_t ievent, Int_t irun, Float_t randNum ) ; // for tests
    void Print(Option_t *option="") const;
    
    void initialize(UInt_t nEvent, UInt_t nRun, Int_t sourceId, UInt_t sequence, Double_t timeStamp);
    
    inline UInt_t getEventId() const { return m_eventId; };
    inline UInt_t getRunId() const { return m_runId; };
	inline Int_t getSourceId() const { return m_sourceId; };
	inline UInt_t getSequence() const { return m_sequence; };

    
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
    /// clear of the array (necessary for converters)
    void clearMcParticleCol()  { m_particleCol->Clear(); };
   
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
    /// clear of the array (necessary for converters)
    void clearMcPositionHitCol()  { m_positionHitCol->Clear(); };
 
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
    /// clear of the array (necessary for converters)
    void clearMcIntegratingHitCol()  { m_integratingHitCol->Clear(); };
   
    /// store a new McTrajectory in the collection
    void addMcTrajectory(McTrajectory *trajectory);
    /// return a MMcTrajectory corresponding to index
    McTrajectory* getMcTrajectory(UInt_t index) const;
    /// return the number of McTrajectories stored in the collection
    inline UInt_t getMcTrajectoryCount() const { 
        return ((m_trajectoryCol) ? m_trajectoryCol->GetEntries() : 0); 
    };
    /// return the full TObjArray containing McPositionHits
    const TObjArray* getMcTrajectoryCol() const { return m_trajectoryCol; };
    /// clear of the array (necessary for converters)
    void clearMcTrajectoryCol()  { m_trajectoryCol->Clear(); };
 
/// time stamp stuff here
    inline Double_t getTimeStamp() const { return m_timeStamp; };

private:
    /// unique event id for this run
    UInt_t m_eventId;
    
    /// Run number
    UInt_t m_runId;

	// identifier of the source
	Int_t m_sourceId;

	/// Sequence number
	UInt_t m_sequence;
    
            /// Time in seconds
    Double_t m_timeStamp; /// this is supposed to coresponed to the following in MCEvent 
        /// Time stamp: use special class to encapsulate type
        ///TimeStamp           m_time;

    /// Collection of McParticles
    /// Must be TObjArray due to TRefArray which grows dynamically
    TObjArray *m_particleCol; //->
    /// static array to allow one-time array creation
    static TObjArray *m_staticParticleCol; //!
    
    /// Collection of McIntegratingHits
    /// must be a TObjArray due to map in McIntegratingHit
    TObjArray *m_integratingHitCol; //->
    /// static array to allow one-time array creation
    static TObjArray *m_staticIntegratingHitCol; //!
    
    /// Collection of McPositionHits
    TObjArray *m_positionHitCol; //->
    /// static array to allow one-time array creation
    static TObjArray *m_staticPositionHitCol; //!
    
    /// Collection of McTrajectories
    TObjArray *m_trajectoryCol; //->
    /// static array to allow one-time array creation
    static TObjArray *m_staticTrajectoryCol; //!
    
    ClassDef(McEvent,6) // Monte Carlo Event Class
};

#endif
