// $Header$
#ifndef GlastEvent_McPositionHit_H
#define GlastEvent_McPositionHit_H 1

// If you wish to introduce the namespace `GlastEvent', uncomment
// the lines commented as `NameSpace'.


// Include files
#include <iostream>
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

class McPositionHit : virtual public ContainedObject {
  public:
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
    /// Packed flags for the internal use.
    unsigned long           m_packedFlags;
};


typedef ObjectVector<McPositionHit> McPositionHitVector;
typedef ObjectList<McPositionHit>   McPositionHitList;

//} // NameSpace GlastEvent


// Inline codes
#include <math.h>
#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/MonteCarlo/McConstants.h"

//namespace GlastEvent { // NameSpace

/// Retrieve cell identifier
inline const VolumeID McPositionHit::volumeID() const
{
  return m_volumeID;
}
/// Update cell identifier
inline void McPositionHit::setVolumeID( VolumeID value )
{
  m_volumeID = value;
}


/// Retrieve entry member
inline const HepPoint3D& McPositionHit::entryPoint() const
{
  return m_entry;
}
inline HepPoint3D& McPositionHit::entryPoint()
{
  return m_entry;
}


/// Update Entry member
inline void McPositionHit::setEntryPoint( const HepPoint3D& value )
{
  m_entry = value;
}


/// Retrieve exit point
inline const HepPoint3D& McPositionHit::exitPoint() const
{
  return m_exit;
}
inline HepPoint3D& McPositionHit::exitPoint()
{
  return m_exit;
}


/// Update exit point
inline void McPositionHit::setExitPoint( const HepPoint3D& value )
{
  m_exit = value;
}


/// Retrieve deposited energy
inline double McPositionHit::depositedEnergy() const
{
  return m_depositedEnergy;
}


/// Update deposited energy
inline void McPositionHit::setDepositedEnergy( double value )
{
  m_depositedEnergy = value;
}

/// Retrieve depositing particle's energy
inline double McPositionHit::particleEnergy() const
{
  return m_particleEnergy;
}

/// Set depositing particle's energy
inline void McPositionHit::setParticleEnergy(double value)
{
  m_particleEnergy = value;
}

/// Retrieve primary-origin flag
inline bool McPositionHit::primaryOrigin() const
{
  using GlastEvent::McConstants::ORIGIN_PRIMARY;
  return m_packedFlags & ORIGIN_PRIMARY;
}
/// Update primary-origin flag
inline void McPositionHit::setPrimaryOrigin( bool value )
{
  using GlastEvent::McConstants::ORIGIN_PRIMARY;
  if (value){
    m_packedFlags |= ORIGIN_PRIMARY;
  } else {
    m_packedFlags &= ~ORIGIN_PRIMARY;
  }
}
/// Retrieve calorimeter-shower-origin flag
inline bool McPositionHit::caloShowerOrigin() const
{
  using GlastEvent::McConstants::ORIGIN_CALOSHOWER;
  return m_packedFlags & ORIGIN_CALOSHOWER;
}
/// Update calorimeter-shower-origin flag
inline void McPositionHit::setCaloShowerOrigin( bool value )
{
  using GlastEvent::McConstants::ORIGIN_CALOSHOWER;
  if (value){
    m_packedFlags |= ORIGIN_CALOSHOWER;
  } else {
    m_packedFlags &= ~ORIGIN_CALOSHOWER;
  }
}

/// Retrieve whether this hit should be digitized
inline bool McPositionHit::needDigi() const
{
  using GlastEvent::McConstants::NEED_DIGI;
  return m_packedFlags & NEED_DIGI;
}
/// Update whether this hit should be digitized
inline void McPositionHit::setNeedDigi( bool value )
{
  using GlastEvent::McConstants::NEED_DIGI;
  if (value){
    m_packedFlags |= NEED_DIGI;
  } else {
    m_packedFlags &= ~NEED_DIGI;
  }
}


/// Retrieve member TOF
inline double McPositionHit::timeOfFlight() const
{
  return m_timeOfFlight;
}


/// Update TOF member
inline void McPositionHit::setTimeOfFlight( double value )
{
  m_timeOfFlight = value;
}


/// Retrieve pointer to McParticle (const or non-const)
inline const McParticle* McPositionHit::mcParticle() const
{
  return m_mcParticle; 
}
inline       McParticle* McPositionHit::mcParticle()
{
  return m_mcParticle; 
}
/// Update pointer to McParticle (by a C++ pointer or a smart reference)
inline void McPositionHit::setMcParticle( McParticle* value )
{
  m_mcParticle = value; 
}
inline void McPositionHit::setMcParticle( SmartRef<McParticle> value )
{ 
  m_mcParticle = value; 
}

//} // NameSpace GlastEvent

#endif    // GlastEvent_McPositionHit_H
