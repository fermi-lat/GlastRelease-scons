#ifndef MCPOSITIONHIT_H
#define MCPOSITIONHIT_H


#include "TObject.h"
#include "TVector.h"
#include "TRef.h"

#include "VolumeIdentifier.h"
#include "McParticle.h"

/* @class McPositionHit
* @brief A ROOT class for storing McPositionHit information.
* This class mirrors the McPositionHit class used in the Gaudi TDS.
* Each McPositionHit contains:
* - VolumeIdentifier
* - Status flags
* - Deposited energy
* - Particle energy
* - Time of Flight
* - Entry and Exit positions
* - Pointer to the McParticle causing the hit
* - Pointer to the origin McParticle
* 
* $Header$
*/

class McPositionHit : public TObject {
    
public:
    
    McPositionHit();
    
    virtual ~McPositionHit();

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    void initialize(Int_t particleId, Double_t edep, 
        const VolumeIdentifier& volId, const TVector3& entry,
        const TVector3& exit, McParticle *mc, McParticle *origin, Double_t pE,
        Double_t tof, UInt_t flags = 0);
    
    const VolumeIdentifier& getVolumeId() const { return m_volumeId; };
    
    const TVector3& getEntryPosition() const { return m_entry; }
    
    const TVector3& getExitPosition() const { return m_exit; };
    
    Int_t getMcParticleId() const { return m_mcParticleId; };

    Double_t getDepositedEnergy() const { return m_depositedEnergy; }
    
    Double_t getParticleEnergy() const { return m_particleEnergy; }
    
    Double_t getTimeOfFlight() const { return m_timeOfFlight; }

    Double_t getDirectionCosine() const;
    
    const McParticle* getMcParticle() const {return (McParticle*)m_mcParticle.GetObject();}
    
    const McParticle* getOriginMcParticle() const { return (McParticle*)m_originMcParticle.GetObject();}
    
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

    TVector3 m_entry;
    TVector3 m_exit;
    
    TRef m_mcParticle;
    TRef m_originMcParticle;
    
    VolumeIdentifier m_volumeId;
    
    ClassDef(McPositionHit,2)  // Monte Carlo PositionHit class
};

#endif