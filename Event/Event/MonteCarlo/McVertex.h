// $Header$
#ifndef GlastEvent_McVertex_H
#define GlastEvent_McVertex_H 1

// If you wish to introduce the namespace `GlastEvent', uncomment
// the lines commented as `NameSpace'.


// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "Gaudi/Kernel/SmartRef.h"
#include "Gaudi/Kernel/SmartRefVector.h"
#include "GlastEvent/TopLevel/Definitions.h"
#include "CLHEP/Geometry/Point3D.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"
// Include all Glast container types here
//   to simplify inlude statements in algorithms
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"



/*!
//------------------------------------------------------------------------------
//
// ClassName:   McVertex
//
// Description: Monte Carlo vertex information
//              (Currently the information copied from the SICB AVMC bank)
//
// The class McVertex uses the Class Library for HEP (CLHEP).
//
// Author:      OZAKI Masanobu
//
// Changes:     M.Frank 04/10/1999 : Proper use of SmartRefs and SmartRefVectors
//              P.Binko 19/10/1999 : Proper accessors of smart references,
//                                   Formating of ASCII output
//              M.Ozaki 2000-12-05 : Modified for GLAST vertex
//                                   Added vertex type and deleted subEvtID
//              M.Ozaki 2001-01-05 : MCVertex -> McVertex
//
//------------------------------------------------------------------------------
 */

//namespace GlastEvent { // NameSpace

// Forward declarations
class McParticle;

class McVertex : virtual public ContainedObject {
  public:
    // vertex type definition
    enum originType {primaryOrigin = 1, daughterOrigin = 2, decayProduct = 3, showerContents = 4, showerBacksplash = 5};

    /// Constructors
    McVertex() :
     m_subEvtID(0),
     m_timeOfFlight(0.),
     m_vertexType(primaryOrigin)
    { }

    /// Destructor (generated)
    virtual ~McVertex() { }


    /// Retrieve initial position
    const HepPoint3D& initialPosition () const;
    /// Retrieve initial position
          HepPoint3D& initialPosition ();
    /// Update initial position
    void setInitialPosition (const HepPoint3D& value);

    /// Retrieve final position
    const HepPoint3D& finalPosition () const;
    /// Retrieve final position
          HepPoint3D& finalPosition ();
    /// Update final position
    void setFinalPosition (const HepPoint3D& value);

    /// retrieve time of flight
    double timeOfFlight () const;
    /// update time of flight
    void setTimeOfFlight (double value);
    /// retrieve vertex type
    originType vertexType () const;
    /// update vertex type
    void setVertexType (originType value);

    /// Retrieve initial 4-momentum
    const HepLorentzVector& initialFourMomentum() const;
    /// Retrieve initial 4-momentum
          HepLorentzVector& initialFourMomentum();
    /// Update initial 4-momentum
    void setInitialFourMomentum( const HepLorentzVector& value );

    /// Retrieve final 4-momentum
    const HepLorentzVector& finalFourMomentum() const;
    /// Retrieve final 4-momentum
          HepLorentzVector& finalFourMomentum();
    /// Update final 4-momentum
    void setFinalFourMomentum( const HepLorentzVector& value );

    /// Retrieve pointer to mother particle (const or non-const)
    const McParticle* motherMcParticle() const;
          McParticle* motherMcParticle();
    /// Update pointer to mother particle (by a C++ pointer or a smart reference)
    void setMotherMcParticle( McParticle* value );
    void setMotherMcParticle( SmartRef<McParticle> value );

    /// Retrieve pointer to vector of daughter particles (const or non-const)
    const SmartRefVector<McParticle>& daughterMcParticles() const;
          SmartRefVector<McParticle>& daughterMcParticles();
    /// Update all daughter particles
    void setDaughterMcParticles( const SmartRefVector<McParticle>& value );
    /// Remove all daughter particles
    void removeDaughterMcParticles();
    /// Add single daughter particle to vector of daughter particles
    ///   (by a C++ pointer or a smart reference)
    void addDaughterMcParticle( McParticle* value );
    void addDaughterMcParticle( SmartRef<McParticle> value );

    /// Retrieve sub event ID
    short subEvtID() const;
    /// Set sub event ID
    void setSubEvtID( short value );

    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;

  private:
    /// Sub-event ID
    short                      m_subEvtID;
    /// Positions:
    /// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/Geometry/HepPoint3D.html">class HepPoint3D</A>
    /// Initial position
    HepPoint3D                 m_initialPosition;
    /// Final position
    HepPoint3D                 m_finalPosition;
    /// Time of fligfht
    double                     m_timeOfFlight;
    /// vertex type
    originType                 m_vertexType;
    /// 4-momentum vectors:
    /// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/Vector/HepLorentzVector.html">class HepLorentzVector</A>
    /// Initial 4-momentum
    HepLorentzVector          m_initialFourMomentum;
    /// Final 4-momentum
    HepLorentzVector          m_finalFourMomentum;
    /// Pointer to mother particle
    SmartRef<McParticle>       m_motherMcParticle;
    /// Vector of pointers to daughter particles
    SmartRefVector<McParticle> m_daughterMcParticles;
};


// Definition of all container types of McVertex
template <class TYPE> class ObjectVector;
typedef ObjectVector<McVertex>     McVertexVector;
template <class TYPE> class ObjectList;
typedef ObjectList<McVertex>       McVertexList;

//} // NameSpace GlastEvent


// Inline codes
#include "GlastEvent/MonteCarlo/McParticle.h"

//namespace GlastEvent { // NameSpace

/// Retrieve initial position
inline const HepPoint3D& McVertex::initialPosition () const
{
  return m_initialPosition;
}
/// Retrieve initial position
inline HepPoint3D& McVertex::initialPosition ()
{
  return m_initialPosition;
}
/// Update initial position
inline void McVertex::setInitialPosition (const HepPoint3D& value)
{
  m_initialPosition = value;
}

/// Retrieve final position
inline const HepPoint3D& McVertex::finalPosition () const
{
  return m_finalPosition;
}
/// Retrieve final position
inline HepPoint3D& McVertex::finalPosition ()
{
  return m_finalPosition;
}
/// Update final position
inline void McVertex::setFinalPosition (const HepPoint3D& value)
{
  m_finalPosition = value;
}

/// retrieve time of flight
inline double McVertex::timeOfFlight () const
{
  return m_timeOfFlight;
}
/// update time of flight
inline void McVertex::setTimeOfFlight (double value)
{
  m_timeOfFlight = value;
}
/// retrieve vertex type
inline McVertex::originType McVertex::vertexType () const
{
  return m_vertexType;
}
/// update vertex type
inline void McVertex::setVertexType (McVertex::originType value)
{
  m_vertexType = value;
}


/// Retrieve pointer to mother particle (const or non-const)
inline const McParticle* McVertex::motherMcParticle() const                             {
  return m_motherMcParticle;
}
inline       McParticle* McVertex::motherMcParticle()
{
  return m_motherMcParticle;
}

/// Update pointer to mother particle (by a C++ pointer or a smart reference)
inline void McVertex::setMotherMcParticle( McParticle* value )
{
  m_motherMcParticle = value;
}
inline void McVertex::setMotherMcParticle( SmartRef<McParticle> value )
{
  m_motherMcParticle = value;
}

/// Retrieve pointer to vector of daughter particles (const or non-const)
inline const SmartRefVector<McParticle>& McVertex::daughterMcParticles() const
{
  return m_daughterMcParticles;
}
inline       SmartRefVector<McParticle>& McVertex::daughterMcParticles()
{
  return m_daughterMcParticles;
}

/// Update all daughter particles
inline void McVertex::setDaughterMcParticles( const SmartRefVector<McParticle>& value )
{
  m_daughterMcParticles = value;
}

/// Remove all daughter particles
inline void McVertex::removeDaughterMcParticles()
{
  m_daughterMcParticles.clear();
}

/// Add single daughter particle to vector of daughter particles
///   (by a C++ pointer or a smart reference)
inline void McVertex::addDaughterMcParticle( McParticle* value )
{
  m_daughterMcParticles.push_back(value);
}
inline void McVertex::addDaughterMcParticle( SmartRef<McParticle> value )
{
  m_daughterMcParticles.push_back(value);
}


/// Retrieve sub event ID
inline short McVertex::subEvtID() const
{
  return m_subEvtID;
}
/// Set sub event ID
inline void McVertex::setSubEvtID( short value )
{
  m_subEvtID = value;
}

//} // NameSpace GlastEvent

#endif    // GlastEvent_McVertex_H
