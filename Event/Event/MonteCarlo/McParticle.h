/** 
* @file McReconAlg.cxx
* @brief Declaration and definition of the TDS object McParticle.
*
*  $Header$
*/
#ifndef GlastEvent_McParticle_H
#define GlastEvent_McParticle_H 1

// If you wish to introduce the namespace `GlastEvent', uncomment
// the lines commented as `NameSpace'.


// Include files
#include <iostream>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRef.h"
#include "GaudiKernel/SmartRefVector.h"

#include "GlastEvent/MonteCarlo/McConstants.h"
#include "GlastEvent/TopLevel/Definitions.h"
#include "GlastEvent/Utilities/ParticleID.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"

// Include all Glast container types here
//   to simplify inlude statements in algorithms
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ObjectList.h"



/*!
//------------------------------------------------------------------------------
//
// ClassName:   McParticle
//  
// Description: The Monte Carlo particle kinematics information
//              (Currently the information copied from the SICB ATMC bank)
//
// The class McParticle uses the Class Library for HEP (CLHEP).
//              
// Author:      OZAKI Masanobu
//
// Changes:     M.Ozaki 2000-12-05 : Based on LHCbEvent's MCParticle rev 1.1.1.2
//              M.Ozaki 2001-01-05 : MCParticle -> McParticle
//
//------------------------------------------------------------------------------
 */


extern const CLID& CLID_McParticle;

namespace mc {  // NameSpace

class McParticle  : virtual public ContainedObject  {
  public:
    typedef int  StdHepId;

    //! status bits modeled on successful SLD scheme
    enum StatusBits{  
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
        Swum =   1<<12  //! this particle was produced by the swimmer
    };



    virtual const CLID& clID() const   { return McParticle::classID(); }
    static const CLID& classID()       { return CLID_McParticle; }
    /// Constructors
    McParticle() :
     m_statusFlags(0)
    {}
    /// Destructor
    virtual ~McParticle() {}

    //! completely initialize a newed object. No other way to set most attributes.
    void init( McParticle* mother, 
        StdHepId id, 
        unsigned int statusBits,
        const HepLorentzVector& initalMomentum,
        const HepLorentzVector& finalMomentum,
        const HepPoint3D& finalPosition);

    /// Retrieve particle identificatio
    ParticleID particleID() const;
    /// Update particle identification
    void setParticleID( ParticleID value );

    /// Retrieve particle property
    StdHepId particleProperty() const;
    /// Update particle identification
    void setParticleProperty( StdHepId value );

    /// Retrieve whether this is a primary particle: true if mother is itself
    bool primaryParticle() const;

    /// Retrieve pointer to the start, end vertex positions 
    const HepPoint3D& initialPosition() const;
    const HepPoint3D& finalPosition() const;

    const HepLorentzVector&  initialFourMomemtum()const;

    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const ;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;

  private:

    /// particle property (such as electron or proton or ....) ID
    StdHepId                  m_particleID;
    /// Bit-field status flag
    unsigned long             m_statusFlags;
    /// Final position
    HepPoint3D                 m_finalPosition;
    /// 4-momentum vectors:
    /// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/Vector/HepLorentzVector.html">class HepLorentzVector</A>
    /// Initial 4-momentum
    HepLorentzVector           m_initialFourMomentum;
    /// Final 4-momentum
    HepLorentzVector           m_finalFourMomentum;
    /// Pointer to mother particle
    SmartRef<McParticle>       m_mother;
    /// Vector of pointers to daughter particles
    SmartRefVector<McParticle> m_daughters;
};


// Definition of all container types of McParticle
//template <class TYPE> class ObjectVector;
//typedef ObjectVector<McParticle>     McParticleVector;

typedef ObjectList<mc::McParticle>       McParticleList;
typedef ObjectList<mc::McParticle>       McParticleCol;

inline StreamBuffer& McParticle::serialize( StreamBuffer& s ) const
{
  ContainedObject::serialize(s);
  return s
#if 0
      << m_particleID
    << m_particleProperty
    << m_subEvtID
    << m_statusFlags
    << m_mcVertex(this)
#endif
    ;
}


/// Serialize the object for reading
inline StreamBuffer& McParticle::serialize( StreamBuffer& s )
{
  ContainedObject::serialize(s);
#if 0
  s
    >> m_particleID
    >> m_particleProperty
    >> m_subEvtID
    >> m_statusFlags
    >> m_mcVertex(this);
#endif
  return s;
}


/// Fill the ASCII output stream
inline std::ostream& McParticle::fillStream( std::ostream& s ) const
{
#if 0
  s << "class McParticle"
    << " (SubEvent:" << m_subEvtID << ")"
    << " :"
    << "\n    Particle ID                = " << m_particleID
    << "\n    Particle Property          = " << m_particleProperty
    << "\n    Sub Event ID               = " << m_subEvtID
    << "\n    McVertex                   = " << m_mcVertex(this);
#endif
  return s;
}


} // NameSpace GlastEvent


#endif    // GlastEvent_McParticle_H
