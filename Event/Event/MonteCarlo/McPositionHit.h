// $Header$
#ifndef GlastEvent_McPositionHit_H
#define GlastEvent_McPositionHit_H 1

// If you wish to introduce the namespace `GlastEvent', uncomment
// the lines commented as `NameSpace'.


// Include files
#include <iostream>
#include <math.h>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "Gaudi/Kernel/SmartRef.h"
#include "GlastEvent/TopLevel/Definitions.h"
#include "CLHEP/Geometry/Point3D.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"
// Include all Gaudi container types here
//   to simplify inlude statements in algorithms
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"
#include "GlastEvent/Utilities/VolumeID.h"


#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/MonteCarlo/McConstants.h"



/*!
//------------------------------------------------------------------------------
//
// ClassName:   McPositionHit
//  
// Description: Monte Carlo class for a position hit such as tracker
//
// The class McPositionHit uses the Class Library for HEP (CLHEP).
//              
// Author:      OZAKI Masanobu
// Changes:     M.Frank 04/10/1999 : Proper use of SmartRefs and SmartRefVectors
//              P.Binko 19/10/1999 : Proper accessors of smart references
//                                   Formating of ASCII output
//              M.Ozaki 2000-12-05 : Modified for GLAST
//              M.Ozaki 2001-01-05 : MCPositionHits -> McPositionHit
//
//------------------------------------------------------------------------------
 */

//namespace GlastEvent { // NameSpace

// Forward declarations
class McParticle;

extern const CLID& CLID_McPositionHit;

class McPositionHit : virtual public ContainedObject {
  public:

    virtual const CLID& clID() const   { return McPositionHit::classID(); }
    static const CLID& classID()       { return CLID_McPositionHit; }
    /// Constructors
    McPositionHit()
      : m_depositedEnergy(0.),
        m_timeOfFlight(0.)
    {}
    /// Destructor
    virtual ~McPositionHit() { }


    /// Retrieve cell identifier
    const VolumeID volumeID() const;
    /// Update cell identifier
    void setVolumeID( VolumeID value );

    /// Retrieve entry member
    const HepPoint3D& entryPoint() const;
          HepPoint3D& entryPoint();
    /// Update Entry member
    void setEntryPoint( const HepPoint3D& value );

    /// Retrieve exit point
    const HepPoint3D& exitPoint() const;
          HepPoint3D& exitPoint();
    /// Update exit point
    void setExitPoint( const HepPoint3D& value );

    /// Retrieve deposited energy
    double depositedEnergy() const;
    /// Update deposited energy
    void setDepositedEnergy( double value );
    /// Retrieve depositing particle's energy
    double particleEnergy() const;
    /// Update depositing particle's energy
    void setParticleEnergy( double value );

    /// Retrieve primary-origin flag
    bool primaryOrigin() const;
    /// Update primary-origin flag
    void setPrimaryOrigin( bool value );
    /// Retrieve calorimeter-shower-origin flag
    bool caloShowerOrigin() const;
    /// Update calorimeter-shower-origin flag
    void setCaloShowerOrigin( bool value );

    /// Retrieve whether this hit should be digitized
    bool needDigi() const;
    /// Update whether this hit should be digitized
    void setNeedDigi( bool value );

    /// Retrieve hit's direction cosine
    double directionCosine() const;

    /// Retrieve member TOF
    double timeOfFlight() const;
    /// Update TOF member
    void setTimeOfFlight( double value );

    /// Retrieve pointer to McParticle (const or non-const)
    const McParticle* mcParticle() const;
          McParticle* mcParticle();
    /// Update pointer to McParticle (by a C++ pointer or a smart reference)
    void setMcParticle( McParticle* value );
    void setMcParticle( SmartRef<McParticle> value );

    /// Retrieve pointer to the ancestor McParticle (const or non-const)
    const McParticle* originMcParticle() const;
          McParticle* originMcParticle();
    /// Update pointer to McParticle (by a C++ pointer or a smart reference)
    void setOriginMcParticle( McParticle* value );
    void setOriginMcParticle( SmartRef<McParticle> value );

    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;

  private:
    /// Volume ID
    VolumeID                m_volumeID;
    /// Entry point
    HepPoint3D              m_entry;
    /// Exit point
    HepPoint3D              m_exit;
    /// Deposited energy
    double                  m_depositedEnergy;
    /// Depositing particle's energy
    double                  m_particleEnergy;
    /// Time of flight
    double                  m_timeOfFlight;
    /// Pointer to McParticle causing the hit
    SmartRef<McParticle>    m_mcParticle;
    /// Pointer to the ancestor McParticle
    SmartRef<McParticle>    m_originMcParticle;
    /// Packed flags for the internal use.
    unsigned long           m_packedFlags;
};


typedef ObjectVector<McPositionHit> McPositionHitVector;
typedef ObjectList<McPositionHit>   McPositionHitList;

/// Serialize the object for writing
inline StreamBuffer& McPositionHit::serialize( StreamBuffer& s ) const
{
  ContainedObject::serialize(s);
  return s
    << m_volumeID
    << m_entry
    << m_exit
    << m_depositedEnergy
    << m_particleEnergy
    << m_timeOfFlight
    << m_mcParticle(this)
    << m_originMcParticle(this)
    << m_packedFlags;
}


/// Serialize the object for reading
inline StreamBuffer& McPositionHit::serialize( StreamBuffer& s )
{
  ContainedObject::serialize(s);
  return s
    >> m_volumeID
    >> m_entry
    >> m_exit
    >> m_depositedEnergy
    >> m_particleEnergy
    >> m_timeOfFlight
    >> m_mcParticle(this)
    >> m_originMcParticle(this)
    >> m_packedFlags;
}


/// Fill the ASCII output stream
inline std::ostream& McPositionHit::fillStream( std::ostream& s ) const
{
  return s
    << "    base class McPositionHit :"
    << "\n        Volume ID             = " << m_volumeID
    << "\n        Entry point (x, y, z) = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_entry.x() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_entry.y() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_entry.z() << " )"
    << "\n        Deposited Energy      = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_depositedEnergy
    << "\n        Particle Energy       = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_particleEnergy
    << "\n        Time of flight        = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_timeOfFlight
    << "\n        Exit point (x, y, z)  = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_exit.x() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_exit.y() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_exit.z() << " )"
    << "\n        McParticle            = " << m_mcParticle(this)
    << "\n        ancestor McParticle   = " << m_originMcParticle(this);
}

#endif    // GlastEvent_McPositionHit_H
