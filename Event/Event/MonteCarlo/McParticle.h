#ifndef Event_McParticle_H
#define Event_McParticle_H 1

// If you wish to introduce the namespace `Event', uncomment
// the lines commented as `NameSpace'.

#include <iostream>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRef.h"
#include "GaudiKernel/SmartRefVector.h"

#include "Event/TopLevel/Definitions.h"
#include "Event/Utilities/ParticleID.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "Event/Utilities/CLHEPStreams.h"

// Include all Glast container types here
//   to simplify inlude statements in algorithms
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ObjectList.h"

/** @class McParticle
 * @brief The Monte Carlo particle kinematics information
 * (Currently the information copied from the SICB ATMC bank)
 *
 *  The class McParticle uses the Class Library for HEP (CLHEP).
 *            
 * @author OZAKI Masanobu
 *
 * Changes:     M.Ozaki 2000-12-05 : Based on LHCbEvent's MCParticle rev 1.1.1.2
 *              M.Ozaki 2001-01-05 : MCParticle -> McParticle
 *
 * $Header$
 */
extern const CLID& CLID_McParticle;

namespace Event {  // NameSpace


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
    //! it will be replaced by the following methods (left here just in the transition)
    void init( McParticle* mother, 
        StdHepId id, 
        unsigned int statusBits,
        const HepLorentzVector& initialMomentum,
        const HepLorentzVector& finalMomentum,
        const HepPoint3D& initialPosition,
        const HepPoint3D& finalPosition,
        const std::string process = "");

    //! Set the initial attributes of the McParticle
    void initialize( McParticle* mother, 
        StdHepId id, 
        unsigned int statusBits,
        const HepLorentzVector& initialMomentum,
        const HepPoint3D& initialPosition,
        const std::string process = "");

    //! Set the final attributes of the McParticle
    void finalize( const HepLorentzVector& finalMomentum,
        const HepPoint3D& finalPosition);


    /// Retrieve particle property
    StdHepId particleProperty() const;

    /// retrieve all of status flags for const object
    unsigned int statusFlags()const;


    /// Retrieve whether this is a primary particle: true if mother is itself
    bool primaryParticle() const;

    /// Retrieve pointer to the start, end vertex positions 
    const HepPoint3D& initialPosition() const;
    const HepPoint3D& finalPosition() const;

    const HepLorentzVector&  initialFourMomentum()const;
    const HepLorentzVector&  finalFourMomentum()const;

    /// access to the mother particle
    const McParticle& mother()const; 

    /// access the process name
    const std::string getProcess()const{return m_process;};


    /// access to the list of daughters
    const SmartRefVector<McParticle>& daughterList()const;

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
    /// Initial position
    HepPoint3D                 m_initialPosition;
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
    /// String with the process name that poduces this particle
    std::string                m_process;
};


// Definition of all container types of McParticle
//template <class TYPE> class ObjectVector;
//typedef ObjectVector<McParticle>     McParticleVector;

typedef ObjectList<McParticle>       McParticleList;
typedef ObjectList<McParticle>       McParticleCol;

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
  return s << "class McParticle"
    << " :"
    << "\n    Particle ID                = " << m_particleID
    << "\n    Status Flags               = " << m_statusFlags
    << "\n    Initial position (x, y, z) = ( "
    << EventFloatFormat (EventFormat::width, EventFormat::precision )
    << m_initialPosition.x() << ","
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_initialPosition.y() << ","
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_initialPosition.z() << " ) "
    << "\n    Final position (x, y, z) = ( "
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_finalPosition.x() << ","
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_finalPosition.y() << ","
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_finalPosition.z() << " ) "
    << "\n    Initial Momentum  (x,y,z,t)   = ( " 
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_initialFourMomentum.x() << ","
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_initialFourMomentum.y() << ","
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_initialFourMomentum.z() << ","
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_initialFourMomentum.t() << " ) "
    << "\n    Final Momentum   (x,y,z,t)    = ( " 
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_finalFourMomentum.x() << ","
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_finalFourMomentum.y() << ","
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_finalFourMomentum.z() << ","
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_finalFourMomentum.t() << " )";
    
    //<< "\n    Mother                     = " << m_mother;
}


} // NameSpace Event


#endif    // Event_McParticle_H
