// $Header$
#ifndef GlastEvent_MCParticle_H
#define GlastEvent_MCParticle_H 1


// Uncomment the following line if you need sub-event ID.
//#define NeedSubEvtID


// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "Gaudi/Kernel/SmartRef.h"
#include "Gaudi/Kernel/SmartRefVector.h"
#include "Gaudi/ParticlePropertySvc/ParticleProperty.h"
#include "GlastEvent/TopLevel/Definitions.h"
#include "GlastEvent/Utilities/ParticleID.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"
// Include all Glast container types here
//   to simplify inlude statements in algorithms
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"


// Forward declarations
class MCVertex;


/*!
//------------------------------------------------------------------------------
//
// ClassName:   MCParticle
//  
// Description: The Monte Carlo particle kinematics information
//              (Currently the information copied from the SICB ATMC bank)
//
// The class MCParticle uses the Class Library for HEP
// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/index.html">CLHEP</A>
//              
// Author:      OZAKI Masanobu
//
// Changes:     M.Ozaki 2000-12-05 : Based on LHCbEvent's MCParticle rev 1.1.1.2
//
//------------------------------------------------------------------------------
 */


class MCParticle  : virtual public ContainedObject  {

public:
    /// Constructors
    MCParticle() :
#ifdef NeedSubEvtID
     m_subEvtID(0),
#endif // NeedSubEvtID
     m_primaryParticle(false),
     m_statusFlags(0)
    {}
    /// Destructor
    virtual ~MCParticle() {}

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
    const MCVertex* originMCVertex() const;
          MCVertex* originMCVertex();
    /// Update pointer to origin vertex (by a C++ pointer or a smart reference)
    void setOriginMCVertex( const MCVertex* value );
    void setOriginMCVertex( const SmartRef<MCVertex> value );

    /// Retrieve pointer to vector of end vertex (const or non-const)
    const MCVertex* endMCVertex() const;
          MCVertex* endMCVertex();
    /// Update pointer to origin vertex (by a C++ pointer or a smart reference)
    void setEndMCVertex( const MCVertex* value );
    void setEndMCVertex( const SmartRef<MCVertex> value );

#ifdef NeedSubEvtID
    /// Retrieve sub event ID 
    short subEvtID() const;
    /// Set sub event ID 
    void setSubEvtID( short value );
#endif // NeedSubEvtID

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
#ifdef NeedSubEvtID
    /// Sub-event ID
    short                     m_subEvtID;
#endif // NeedSubEvtID
    /// Bit-field status flag
    unsigned long             m_statusFlags;
    /// Pointer to origin vertex
    SmartRef<MCVertex>        m_originMCVertex;
    /// Vector of pointers to end vertex
    SmartRef<MCVertex>        m_endMCVertex;

    /// Constant(s) for internal use.  Should be moved to .cpp file when
    /// the namespace "GlastEvent" is introduced.
    const unsigned long       STATUS_PRIMARY = 1<<0;
};


// Definition of all container types of MCParticle
template <class TYPE> class ObjectVector;
typedef ObjectVector<MCParticle>     MCParticleVector;
template <class TYPE> class ObjectList;
typedef ObjectList<MCParticle>       MCParticleList;


// Inline codes
#include "GlastEvent/MonteCarlo/MCParticle.cpp"


#endif    // GlastEvent_MCParticle_H
