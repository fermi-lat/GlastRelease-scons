#ifndef MCPARTICLE_H
#define MCPARTICLE_H

#include "TObject.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "TRef.h"
#include "TRefArray.h"

/**
* @class McParticle
* @brief GLAST Monte Carlo particle class.  Designed originally for 2001 balloon flight.
*  Inheritance information given for parentage only.
*
*  \author Daniel Flath and Heather Kelly
*  
*/
class McParticle: public TObject {
    
public:
    //! status bits modeled on successful SLD scheme
    typedef enum StatusBits{  
        DECAYED =1 ,    //! Decayed by generator
            DECAYFLT=1<<1,  //! Decayed in flight by swimmer
            MISSED=  1<<2,  //! Does not hit detector
            NOINTER =1<<3,  //! Traverses detector w/o interacting 
            STOPPED =1<<4,  //! Energy below cut; other bits may say why 
            INTERACT=1<<5,  //! Interacted, no further decision to be made
            INTSHDEP=1<<6,  //! Interacted, further decision depends on ! selection of shower deposition  
            PRIMARY =1<<7,  //! primary particle 
            SWERROR =1<<8,  //! Error occurred in swimming the track 
            SW2MNYST=1<<9,  //! Swim aborted: too many steps (ISTOP=99)  
            WOUTOFT =1<<10, //! Swim aborted: out of sensitive time of ! detector (ISTOP=4) 
            NOTTRACK=1<<11, //! Not tracked by user request 
            Swum =   1<<12,  //! this particle was produced by the swimmer
            CALOSHOWER = 1 << 13,
            ORIGIN_PRIMARY = 1 << 14,
            ORIGIN_CALOSHOWER  = 1 << 15,
            NEED_DIGI  = 1 << 16
    };
    
    McParticle();

    McParticle(const McParticle& p);
    
    virtual ~McParticle();

    /// clear lists, free pointers, etc., after read from / write to file
    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    //! completely initialize a newed object. No other way to set most attributes.
    void initialize( McParticle* mother, 
        Int_t id, 
        UInt_t statusBits,
        const TLorentzVector& initalMomentum,
        const TLorentzVector& finalMomentum,
        const TVector3& finalPosition);
    
    Int_t getParticleId() const { return m_particleId;};    
    
    /// return an McParticle pointer to the mother particle
    const McParticle* getMother();

    /// return a daughter McParticle corresponding to index
    const McParticle* getDaughter(Int_t index) const;

    /// return the full list of daughters
    const TRefArray& getDaughterList() const { return m_daughters; };

    Int_t getParticleProperty() const;
    
    UInt_t getStatusFlags() const { return m_statusFlags; };

    /// Retrieve whether this is a primary particle
    Bool_t primaryParticle() const;
    
    /// Retrieve pointer to the start vertex positions
    const TVector3& getInitialPosition() const;
    
    /// Retrieve pointer to end vertex position
    const TVector3& getFinalPosition() const;
    
    const TLorentzVector&  getInitialFourMomentum() const;

    const TLorentzVector& getFinalFourMomentum() const;
    
    
private:
    
    /// particle property (such as electron or proton or ....) ID
    Int_t m_particleId;
    /// Bit-field status flag
    UInt_t m_statusFlags;
    /// Final position
    TVector3 m_finalPosition;
    /// Initial 4-momentum
    TLorentzVector m_initialFourMomentum;
    /// Final 4-momentum
    TLorentzVector m_finalFourMomentum;
    /// Pointer to mother particle
    TRef m_mother;
    /// Reference array of daughter particles
    TRefArray m_daughters;

    ClassDef(McParticle,2) // Monte Carlo Particle Class
};

#endif
