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
#include "CLHEP/Vector/LorentzVector.h"
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

extern const CLID& CLID_McVertex;

class McVertex : virtual public ContainedObject {
  public:

    virtual const CLID& clID() const   { return McVertex::classID(); }
    static const CLID& classID()       { return CLID_McVertex; }
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

    /// Retrieve pointer to the pair particle (const or non-const)
    const McParticle* mcParticle() const;
          McParticle* mcParticle();
    /// Update pointer to the pair particle (by a C++ pointer or a smart reference)
    void setMcParticle( McParticle* value );
    void setMcParticle( SmartRef<McParticle> value );

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
    HepLorentzVector           m_initialFourMomentum;
    /// Final 4-momentum
    HepLorentzVector           m_finalFourMomentum;
    /// The pair McParticle
    SmartRef<McParticle>       m_mcParticle;
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


//} // NameSpace GlastEvent

#endif    // GlastEvent_McVertex_H
