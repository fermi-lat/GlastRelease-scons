// $Header$
#ifndef GlastEvent_McParticle_H
#define GlastEvent_McParticle_H 1

// If you wish to introduce the namespace `GlastEvent', uncomment
// the lines commented as `NameSpace'.


// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "Gaudi/Kernel/SmartRef.h"
#include "Gaudi/Kernel/SmartRefVector.h"
#include "GlastEvent/MonteCarlo/McVertex.h"
#include "GlastEvent/MonteCarlo/McConstants.h"
#if 0 // *** FIXME!! ***
    // Gaudi v8's ParticleProperty class seems to have some problems
    // around the message stream, so temporarily commented-out.
#include "Gaudi/ParticlePropertySvc/ParticleProperty.h"
#else // 0
    /// An ad-hoc workaround...
#include "Gaudi/Kernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/Definitions.h"
    struct ParticleProperty {
        int m_dummy;
        typedef ParticleProperty PP;
        friend StreamBuffer& operator<< ( StreamBuffer& s, const PP& obj ){
            return s << obj.m_dummy;
        }
        /// Serialize the object for reading
        friend StreamBuffer& operator>> ( StreamBuffer& s, PP& obj ){
            return s >> obj.m_dummy;
        }
        /// Output operator (ASCII)
        friend std::ostream& operator<< ( std::ostream& s, const PP& obj ){
            return obj.fillStream(s);
        }
        /// Fill the output stream (ASCII)
        std::ostream& fillStream( std::ostream& s ) const{
            return s << "class ParticleProperty : "
              << GlastEventField( GlastEvent::field4 )
              << m_dummy;
        }
    };
#endif // 0
#include "GlastEvent/TopLevel/Definitions.h"
#include "GlastEvent/Utilities/ParticleID.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"
// Include all Glast container types here
//   to simplify inlude statements in algorithms
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"



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

//namespace GlastEvent {  // NameSpace

// Forward declarations
class McVertex;

class McParticle  : virtual public ContainedObject  {
  public:
    /// Constructors
    McParticle() :
     m_subEvtID(0),
//   m_primaryParticle(false),
     m_statusFlags(0)
    {}
    /// Destructor
    virtual ~McParticle() {}

    /// Retrieve particle identification
    ParticleID particleID() const;
    /// Update particle identification
    void setParticleID( ParticleID value );

    /// Retrieve particle property
    ParticleProperty particleProperty() const;
    /// Update particle identification
    void setParticleProperty( ParticleProperty value );

    /// Retrieve whether this is a primary particle
    bool primaryParticle() const;
    /// Set whether this is a primary particle
    void setPrimaryParticleFlag( bool value );

    /// Retrieve pointer to the vertex (const or non-const)
    const McVertex* mcVertex() const;
          McVertex* mcVertex();
    /// Update pointer to origin vertex (by a C++ pointer or a smart reference)
    void setMcVertex( McVertex* value );
    void setMcVertex( SmartRef<McVertex> value );

    /// Retrieve sub event ID 
    short subEvtID() const;
    /// Set sub event ID 
    void setSubEvtID( short value );

    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const ;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;

  private:
    /// Particle ID
    ParticleID                m_particleID;
    /// We need particle property (such as electron or proton or ....)
    ParticleProperty          m_particleProperty;
    /// Sub-event ID
    short                     m_subEvtID;
    /// Bit-field status flag
    unsigned long             m_statusFlags;
    /// Pointer to the McVertex
    SmartRef<McVertex>        m_mcVertex;
};


// Definition of all container types of McParticle
template <class TYPE> class ObjectVector;
typedef ObjectVector<McParticle>     McParticleVector;
template <class TYPE> class ObjectList;
typedef ObjectList<McParticle>       McParticleList;

//} // NameSpace GlastEvent


#endif    // GlastEvent_McParticle_H
