// $Header$
#ifndef GlastEvent_MCVertex_H
#define GlastEvent_MCVertex_H 1


// Uncomment the following line if you need sub-event ID.
//#define NeedSubEvtID


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


// Forward declarations
class MCParticle;  // aka MCParticleKinematics


/*!
//------------------------------------------------------------------------------
//
// ClassName:   MCVertex
//
// Description: Monte Carlo vertex information
//              (Currently the information copied from the SICB AVMC bank)
//
// The class MCVertex uses the Class Library for HEP
// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/index.html">CLHEP</A>
//
// Author:      OZAKI Masanobu
//
// Changes:     M.Frank 04/10/1999 : Proper use of SmartRefs and SmartRefVectors
//              P.Binko 19/10/1999 : Proper accessors of smart references,
//                                   Formating of ASCII output
//              M.Ozaki 2000-12-05 : Modified for GLAST vertex
//                                   Added vertex type and deleted subEvtID
//
//------------------------------------------------------------------------------
 */


class MCVertex : virtual public ContainedObject {

public:
    // vertex type definition
    enum originType {primaryOrigin = 1, daughterOrigin = 2, decayProduct = 3, showerContents = 4, showerBacksplash = 5};

    /// Constructors
    MCVertex() :
#ifdef NeedSubEvtID
     m_subEvtID(0),
#endif // NeedSubEvtID
     m_timeOfFlight(0.),
     m_vertexType(primaryOrigin)
    { }

    /// Destructor (generated)
    virtual ~MCVertex() { }


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
    const MCParticle* motherMCParticle() const;
          MCParticle* motherMCParticle();
    /// Update pointer to mother particle (by a C++ pointer or a smart reference)
    void setMotherMCParticle( const MCParticle* value );
    void setMotherMCParticle( const SmartRef<MCParticle> value );

    /// Retrieve pointer to vector of daughter particles (const or non-const)
    const SmartRefVector<MCParticle>& daughterMCParticles() const;
          SmartRefVector<MCParticle>& daughterMCParticles();
    /// Update all daughter particles
    void setDaughterMCParticles( const SmartRefVector<MCParticle>& value );
    /// Remove all daughter particles
    void removeDaughterMCParticles();
    /// Add single daughter particle to vector of daughter particles
    ///   (by a C++ pointer or a smart reference)
    void addDaughterMCParticle( const MCParticle* value );
    void addDaughterMCParticle( const SmartRef<MCParticle> value );

#ifdef NeedSubEvtID
    /// Retrieve sub event ID
    short subEvtID() const;
    /// Set sub event ID
    void setSubEvtID( short value );
#endif // NeedSubEvtID

    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;

private:
#ifdef NeedSubEvtID
    /// Sub-event ID
    short                      m_subEvtID;
#endif // NeedSubEvtID
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
    SmartRef<MCParticle>       m_motherMCParticle;
    /// Vector of pointers to daughter particles
    SmartRefVector<MCParticle> m_daughterMCParticles;
};


// Definition of all container types of MCVertex
template <class TYPE> class ObjectVector;
typedef ObjectVector<MCVertex>     MCVertexVector;
template <class TYPE> class ObjectList;
typedef ObjectList<MCVertex>       MCVertexList;


// Inline codes
#include "GlastEvent/MonteCarlo/MCVertex.cpp"


#endif    // GlastEvent_MCVertex_H
