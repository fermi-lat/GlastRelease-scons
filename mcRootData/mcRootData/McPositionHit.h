#ifndef MCPOSITIONHIT_H
#define MCPOSITIONHIT_H


#include "TObject.h"
#include "TVector.h"
#include "TRef.h"

#include "VolumeIdentifier.h"
#include "McParticle.h"

/* @class McPositionHit
* @brief A ROOT class for storing McPositionHit information.
*
* This class mirrors the McPositionHit class used in the Gaudi TDS.
* Each McPositionHit contains:
* - VolumeIdentifier
* - Status flags
* - Deposited energy
* - Particle energy
* - Time of Flight
* - Entry and Exit positions, in both local and global coordinates
* - Id of McParticle
* - Optional Pointer to McParticle
* - Id of Origin McParticle
* - Optional Pointer to the origin McParticle
* 
* $Header$
*/

class McPositionHit : public TObject {
    
public:
    
    McPositionHit();
    
    virtual ~McPositionHit();

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    void initialize(Int_t particleId, 
        Double_t edep, const VolumeIdentifier& volId, 
        const TVector3& entry, const TVector3& exit, 
        McParticle *mc, McParticle *origin, Double_t pE,
        Double_t tof, UInt_t flags = 0);

    void initialize(Int_t mcParticleId, Int_t originParticleId, 
        Double_t edep, const VolumeIdentifier& volId, 
        const TVector3& entry, const TVector3& exit,
        const TVector3& gEntry, const TVector3& gExit, 
        McParticle *mc, McParticle *origin, Double_t pE,
        Double_t tof, UInt_t flags = 0);
    
    const VolumeIdentifier& getVolumeId() const { return m_volumeId; };
    
    const TVector3& getEntryPosition() const { return m_entry; }

    const TVector3& getGlobalEntryPosition() const { return m_globalEntry; }
    
    const TVector3& getExitPosition() const { return m_exit; };

    const TVector3& getGlobalExitPosition() const { return m_globalExit; }
    
    Int_t getMcParticleId() const { return m_mcParticleId; };

    Int_t getOriginMcParticleId() const { return m_originMcParticleId; };

    Double_t getDepositedEnergy() const { return m_depositedEnergy; }
    
    Double_t getParticleEnergy() const { return m_particleEnergy; }
    
    Double_t getTimeOfFlight() const { return m_timeOfFlight; }

    Double_t getDirectionCosine() const;
        
    const McParticle* getMcParticle() const { return (McParticle*)m_mcParticle.GetObject();};

    const McParticle* getOriginMcParticle() const { return (McParticle*)m_originMcParticle.GetObject(); };
    
    /// Retrieve whether this hit should be digitized
    Bool_t needDigi() const;

private:

    Double_t m_depositedEnergy;
    Double_t m_particleEnergy;
    Double_t m_timeOfFlight;
   
    /// Packed flags for the internal use.
    UInt_t m_statusFlags;
    /// ID of the McParticle causing the hit
    Int_t m_mcParticleId;
    /// ID of the origin McParticle
    Int_t m_originMcParticleId;

    /// local coordinates
    TVector3 m_entry;
    /// local coordinates
    TVector3 m_exit;

    // global coordinates
    TVector3 m_globalEntry;
    // global coordinates
    TVector3 m_globalExit;
    
    TRef m_mcParticle;
    TRef m_originMcParticle;
    
    VolumeIdentifier m_volumeId;
    
    ClassDef(McPositionHit,4)  // Monte Carlo PositionHit class
};

#endif