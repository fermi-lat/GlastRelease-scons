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

    /// Retrieve pointer to origin vertex (const or non-const)
    const McVertex* originMcVertex() const;
          McVertex* originMcVertex();
    /// Update pointer to origin vertex (by a C++ pointer or a smart reference)
    void setOriginMcVertex( McVertex* value );
    void setOriginMcVertex( SmartRef<McVertex> value );

    /// Retrieve pointer to vector of end vertex (const or non-const)
    const McVertex* endMcVertex() const;
          McVertex* endMcVertex();
    /// Update pointer to origin vertex (by a C++ pointer or a smart reference)
    void setEndMcVertex( McVertex* value );
    void setEndMcVertex( SmartRef<McVertex> value );

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
    /// Pointer to origin vertex
    SmartRef<McVertex>        m_originMcVertex;
    /// Pointer to end vertex
    SmartRef<McVertex>        m_endMcVertex;
};


// Definition of all container types of McParticle
template <class TYPE> class ObjectVector;
typedef ObjectVector<McParticle>     McParticleVector;
template <class TYPE> class ObjectList;
typedef ObjectList<McParticle>       McParticleList;

//} // NameSpace GlastEvent



// Inline codes
#include "GlastEvent/MonteCarlo/McVertex.h"
#include "GlastEvent/MonteCarlo/McConstants.h"

//namespace GlastEvent { // NameSpace

/// Retrieve particle identification
inline ParticleID McParticle::particleID() const
{
  return m_particleID;
}


/// Update particle identification
inline void McParticle::setParticleID( ParticleID value )
{
  m_particleID = value;
}


/// Retrieve particle property
inline ParticleProperty McParticle::particleProperty() const
{
  return m_particleProperty;
}


/// Update particle property
inline void McParticle::setParticleProperty( ParticleProperty value )
{
  m_particleProperty = value;
}


/// Retrieve whether this is a primary particle
inline bool McParticle::primaryParticle() const
{
  using GlastEvent::McConstants::PRIMARY;
  return m_statusFlags & PRIMARY;
}


/// Set whether this is a primary particle
inline void McParticle::setPrimaryParticleFlag( bool value )
{
  using GlastEvent::McConstants::PRIMARY;
  if (value){
    m_statusFlags |= PRIMARY;
  } else {
    m_statusFlags &= ~PRIMARY;
  }
}


/// Retrieve pointer to origin vertex (const or non-const)
inline const McVertex* McParticle::originMcVertex() const
{ 
//  McVertex* ret = m_originMcVertex;
//  return ret;
  return m_originMcVertex; 
}
inline       McVertex* McParticle::originMcVertex()
{ 
  return m_originMcVertex;
}


/// Update pointer to origin vertex (by a C++ pointer or a smart reference)
inline void McParticle::setOriginMcVertex( McVertex* value )
{ 
  m_originMcVertex = value; 
}
inline void McParticle::setOriginMcVertex( SmartRef<McVertex> value )
{ 
  m_originMcVertex = value; 
}


/// Retrieve pointer to vector of end vertex (const or non-const)
inline const McVertex* McParticle::endMcVertex() const
{ 
  return m_endMcVertex; 
}
inline       McVertex* McParticle::endMcVertex()
{ 
  return m_endMcVertex; 
}


/// Update pointer to end vertex (by a C++ pointer or a smart reference)
inline void McParticle::setEndMcVertex( McVertex* value )
{ 
  m_endMcVertex = value; 
}
inline void McParticle::setEndMcVertex( SmartRef<McVertex> value )
{ 
  m_endMcVertex = value; 
}


/// Retrieve sub event ID 
inline short McParticle::subEvtID() const
{
  return m_subEvtID;
}


/// Set sub event ID 
inline void McParticle::setSubEvtID( short value )
{
  m_subEvtID = value;
}

//} // NameSpace GlastEvent

#endif    // GlastEvent_McParticle_H
