// $Header$
#ifndef GlastEvent_McPositionHit_H
#define GlastEvent_McPositionHit_H 1


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


// Forward declarations
class MCParticle;


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

    /// Retrieve pointer to MCParticle (const or non-const)
    const MCParticle* mcParticle() const;
          MCParticle* mcParticle();
    /// Update pointer to MCParticle (by a C++ pointer or a smart reference)
    void setMCParticle( const MCParticle* value );
    void setMCParticle( const SmartRef<MCParticle> value );

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
    /// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/Geometry/HepPoint3D.html">class HepPoint3D</A>
    HepPoint3D              m_entry;
    /// Exit point
    HepPoint3D              m_exit;
    /// Deposited energy
    double                  m_depositedEnergy;
    /// Depositing particle's energy
    double                  m_particleEnergy;
    /// Time of flight
    double                  m_timeOfFlight;
    /// Pointer to MCParticle causing the hit
    SmartRef<MCParticle>    m_mcParticle;
    /// Packed flags for the internal use.
    unsigned long           m_packedFlags;

    /// Constants for packedFlags.  Should be moved to .cpp file
    /// when the namespace "GlastEvent" is introduced.
    const unsigned long     ORIGIN_PRIMARY     = 1 << 0;
    const unsigned long     ORIGIN_CALOSHOWER  = 1 << 1;
    const unsigned long     NEED_DIGI          = 1 << 8;
};


typedef ObjectVector<McPositionHit> McPositionHitVector;
typedef ObjectList<McPositionHit>   McPositionHitList;


// Inline codes
#include "GlastEvent/MonteCarlo/McPositionHit.cpp"


#endif    // GlastEvent_McPositionHit_H
