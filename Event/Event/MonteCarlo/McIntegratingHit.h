// $Header$
#ifndef GlastEvent_McIntegratingHit_H
#define GlastEvent_McIntegratingHit_H 1

// If you wish to introduce the namespace `GlastEvent', uncomment
// the lines commented as `NameSpace'.


// Include files
#include <iostream>
#include <map>
#include "CLHEP/Geometry/Point3D.h"
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "Gaudi/Kernel/SmartRefVector.h"
#include "GlastEvent/TopLevel/Definitions.h"
#include "GlastEvent/Utilities/VolumeID.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"
// Include all Glast container types here
//   to simplify inlude statements in algorithms
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"



/*!
//------------------------------------------------------------------------------
//
// ClassName:   McIntegratingHit
//  
// Description: Monte Carlo integrating hits such as calorimeter
//
// Author:      OZAKI Masanobu
//
// Changes:     M.Frank 04/10/1999 : Proper use of SmartRefs and SmartRefVectors
//              P.Binko 19/10/1999 : Proper accessors of smart references,
//                                   Formating of ASCII output
//              M.Ozaki 2000-12-07 : Modified for GLAST
//              M.Ozaki 2001-01-05 : MCIntegratingHits -> McIntegratingHit
//
//------------------------------------------------------------------------------
 */

//namespace GlastEvent {  // NameSpace

// Forward declarations
class McParticle;

class McIntegratingHit : virtual public ContainedObject {
  public:
    /// McParticle -> deposited energy map
    typedef std::map<SmartRef<McParticle>,double>  energyDepositMap;

    /// Constructors
    McIntegratingHit() : m_packedFlags(0)
    {}
    /// Destructor
    ~McIntegratingHit(){}


    /// Retrieve cell identifier
    const VolumeID volumeID() const;
    /// Update cell identifier
    void setVolumeID( VolumeID value );

    /// Retrieve energy
    double totalEnergy() const;
    /// Retrieve the energy-weighted first moments of the position
    const HepPoint3D moment1 () const;
          HepPoint3D moment1 ();
    /// Retrieve the energy-weighted second moments of the position
    const HepPoint3D moment2 () const;
          HepPoint3D moment2 ();

    /// Retrieve itemized energy
    const energyDepositMap& itemizedEnergy() const;
          energyDepositMap& itemizedEnergy();
    /// Update all energyInfos
    void setEnergyItems( const energyDepositMap& value );
    /// Remove all energyInfos
    void clearEnergyItems();
    /// Add single energyInfo to energyDepositMap
    void addEnergyItem( const double& energy, McParticle* t, const HepPoint3D& position );
    void addEnergyItem( const double& energy, SmartRef<McParticle> t, const HepPoint3D& position );

    /// Retrieve primary-origin flag
    bool primaryOrigin() const;
    /// Update primary-origin flag
    void setPrimaryOrigin( bool value );

    /// Retrieve whether this hit should be digitized
    bool needDigi() const;
    /// Update whether this hit should be digitized
    void setNeedDigi( bool value );

    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;

  private:
    /// Cell identifier
    VolumeID                      m_volumeID;
    /// Vector of Energy information that consists of deposited energy and the mother McParticle
    energyDepositMap              m_energyItem;
    /// total deposited energy: set automatically when m_energyInfo is modified.
    double                        m_totalEnergy;
    /// Energy-weighted_first_moments_of_the_position * number_of_energy_deposition
    HepPoint3D                    m_moment1seed;
    /// Energy-weighted_second_moments_of_the_position * number_of_energy_deposition
    HepPoint3D                    m_moment2seed;
    /// Packed flags for particle property
    unsigned long                 m_packedFlags;
};

  
// Definition of all container types of McIntegratingHit
template <class TYPE> class ObjectVector;
typedef ObjectVector<McIntegratingHit>     McIntegratingHitVector;
template <class TYPE> class ObjectList;
typedef ObjectList<McIntegratingHit>       McIntegratingHitList;

//} // NameSpace GlastEvent


// Inline codes
#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/MonteCarlo/McVertex.h"
#include "GlastEvent/MonteCarlo/McConstants.h"

//namespace GlastEvent {  // NameSpace

/// Retrieve volume identifier
inline const VolumeID McIntegratingHit::volumeID() const
{
  return m_volumeID;
}


/// Update volume identifier
inline void McIntegratingHit::setVolumeID( VolumeID value )
{
  m_volumeID = value;
}


/// Retrieve energy
inline double McIntegratingHit::totalEnergy() const
{
  return m_totalEnergy;
}


/// Retrieve the energy-weighted first moments of the position
inline const HepPoint3D McIntegratingHit::moment1 () const
{
    return m_moment1seed * (1./m_totalEnergy);
}
inline HepPoint3D McIntegratingHit::moment1 ()
{
    return m_moment1seed * (1./m_totalEnergy);
}


/// Retrieve the energy-weighted second moments of the position
inline const HepPoint3D McIntegratingHit::moment2 () const
{
    return m_moment2seed * (1./m_totalEnergy);
}
inline HepPoint3D McIntegratingHit::moment2 ()
{
    return m_moment2seed * (1./m_totalEnergy);
}


/// Retrieve itemized energy
inline const McIntegratingHit::energyDepositMap& McIntegratingHit::itemizedEnergy() const
{
  return m_energyItem;
}

inline       McIntegratingHit::energyDepositMap& McIntegratingHit::itemizedEnergy()
{
  return m_energyItem;
}


/// Add an energyItem
inline void McIntegratingHit::addEnergyItem(const double& energy, McParticle* t, const HepPoint3D& position)
{
    m_energyItem[t] += energy;

    HepPoint3D        position2 = HepPoint3D(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
    m_totalEnergy      += energy;
    m_moment1seed += energy * position;
    m_moment2seed += energy * position2;
}


/// Add an energyItem
inline void McIntegratingHit::addEnergyItem(const double& energy, SmartRef<McParticle> t, const HepPoint3D& position)
{
    m_energyItem[t] += energy;

    HepPoint3D        position2 = HepPoint3D(position.x()*position.x(), position.y()*position.y(), position.z()*position.z());
    m_totalEnergy      += energy;
    m_moment1seed += energy * position;
    m_moment2seed += energy * position2;
}


/// Retrieve primary-origin flag
inline bool McIntegratingHit::primaryOrigin() const
{
    using GlastEvent::McConstants::ORIGIN_PRIMARY;
    return m_packedFlags & ORIGIN_PRIMARY;
}


/// Update primary-origin flag
inline void McIntegratingHit::setPrimaryOrigin( bool value )
{
    using GlastEvent::McConstants::ORIGIN_PRIMARY;
    if (value){
        m_packedFlags |= ORIGIN_PRIMARY;
    } else {
        m_packedFlags &= ~ORIGIN_PRIMARY;
    }
}


/// Retrieve whether this hit should be digitized
inline bool McIntegratingHit::needDigi() const
{
    using GlastEvent::McConstants::NEED_DIGI;
    return m_packedFlags & NEED_DIGI;
}


/// Update whether this hit should be digitized
inline void McIntegratingHit::setNeedDigi( bool value )
{
    using GlastEvent::McConstants::NEED_DIGI;
    if (value){
        m_packedFlags |= NEED_DIGI;
    } else {
        m_packedFlags &= ~NEED_DIGI;
    }
}

//} // NameSpace GlastEvent

#endif // GlastEvent_McIntegratingHit_H
