// $Header$

#ifndef _H_GlastHitsEvt_MCTrack_v0_
#define _H_GlastHitsEvt_MCTrack_v0_ 1


// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "Gaudi/Kernel/SmartRef.h"
#include "Gaudi/Kernel/SmartRefVector.h"
#include "GlastEvent/TopLevel/Definitions.h"
#include "GlastEvent/Utilities/ParticleID.h"
#include "CLHEP/Vector/LorentzVector.h"
// Ian#include "Gaudi/Kernel/StreamBuffer.h"

#include "GlastEvent/Utilities/CLHEPStreams.h"
// Include all LHCb container types here
//   to simplify inlude statements in algorithms
//#include "LHCbEvent/TopLevel/ObjectVector.h"
//#include "LHCbEvent/TopLevel/ObjectList.h"


// Forward declarations
//class MCVertex;


// Externals 
extern const CLID& CLID_MCTrack;


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
// Author:      Pavel Binko
//
// Changes:     M.Frank 04/10/1999 : Proper use of SmartRefs and SmartRefVectors
//              P.Binko 19/10/1999 : Proper accessors of smart references
//                                   Formating of ASCII output
//
//------------------------------------------------------------------------------


class MCTrack  : virtual public ContainedObject  {

public:
  /// Constructors
//Ian  MCParticle() : m_oscillationFlag(false)     { } 
  /// Destructor
  virtual ~MCTrack()                       { }

  /// Retrieve pointer to class defininition structure
  virtual const CLID& clID() const            { return MCTrack::classID(); }
  static const CLID& classID()                { return CLID_MCTrack; }

  /// Retrieve virtual mass
  double virtualMass()                        { return m_fourMomentum.m(); }

  /// Retrieve 4-momentum vector
  HepLorentzVector& fourMomentum()                      { return m_fourMomentum; }
  /// Retrieve 4-momentum vector
  const HepLorentzVector& fourMomentum() const          { return m_fourMomentum; }
  /// Update 4-momentum vector  
  void setFourMomentum( const HepLorentzVector& value ) { m_fourMomentum = value; }

  /// Retrieve particle identification
  ParticleID particleID() const               { return m_particleID; }
  /// Update particle identification
  void setParticleID( ParticleID value )      { m_particleID = value; }


/* MCVertex needs to be written before we can include it in MCTrack

  /// Retrieve pointer to origin vertex (const or non-const)
  const MCVertex* originMCVertex() const;
        MCVertex* originMCVertex();
		/// Update pointer to origin vertex (by a C++ pointer or a smart reference)
  void setOriginMCVertex( MCVertex* value );
  void setOriginMCVertex( SmartRef<MCVertex> value );

  /// Retrieve pointer to vector of decay vertices (const or non-const)
  const SmartRefVector<MCVertex>& decayMCVertices() const;
        SmartRefVector<MCVertex>& decayMCVertices();
  /// Update all decay vertices
  void setDecayMCVertices( const SmartRefVector<MCVertex>& value );
  /// Remove all decay vertices
  void removeDecayMCVertices();
  /// Add single decay vertex to vector of decay vertices
  ///   (by a C++ pointer or a smart reference)
  void addDecayMCVertex( MCVertex* value );
  void addDecayMCVertex( SmartRef<MCVertex> value );

  /// Retrieve oscillation Flag
  bool oscillationFlag() const                { return m_oscillationFlag; }
  /// Set oscillation flag
  void setOscillationFlag( bool value )       { m_oscillationFlag = value; }
*/
  /// Serialize the object for writing
  virtual StreamBuffer& serialize( StreamBuffer& s ) const ;
  /// Serialize the object for reading
  virtual StreamBuffer& serialize( StreamBuffer& s );
  /// Fill the ASCII output stream
  virtual std::ostream& fillStream( std::ostream& s ) const;
  

private:
  /// 4-momentum vector
  /// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/Vector/HepLorentzVector.html">class HepLorentzVector</A>
  HepLorentzVector          m_fourMomentum;
  /// Particle ID
  ParticleID                m_particleID;

/* Need MCVertex
  /// Pointer to origin vertex
  SmartRef<MCVertex>        m_originMCVertex;
  /// Vector of pointers to decay vertices
  SmartRefVector<MCVertex>  m_decayMCVertices;
  */

};


//
// Inline code must be outside the class definition
//
//#include "LHCbEvent/MonteCarlo/MCVertex.h"


/// Retrieve pointer to origin vertex (const or non-const)

/* More code that requires MCVertex

inline const MCVertex* MCParticle::originMCVertex() const             { 
  return m_originMCVertex; 
}
inline       MCVertex* MCParticle::originMCVertex()                   { 
  return m_originMCVertex;
}
// Update pointer to origin vertex (by a C++ pointer or a smart reference)
inline void MCParticle::setOriginMCVertex( MCVertex* value )          { 
  m_originMCVertex = value; 
}
inline void MCParticle::setOriginMCVertex( SmartRef<MCVertex> value ) { 
  m_originMCVertex = value; 
}


/// Retrieve pointer to vector of decay vertices (const or non-const)
inline const SmartRefVector<MCVertex>& MCParticle::decayMCVertices() const          { 
  return m_decayMCVertices; 
}
inline       SmartRefVector<MCVertex>& MCParticle::decayMCVertices()                { 
  return m_decayMCVertices; 
}
/// Update all decay vertices
inline void MCParticle::setDecayMCVertices( const SmartRefVector<MCVertex>& value ) { 
  m_decayMCVertices = value; 
}
/// Remove all decay vertices
inline void MCParticle::removeDecayMCVertices()                                     { 
  m_decayMCVertices.clear(); 
}
/// Add single decay vertex to vector of decay vertices
///   (by a C++ pointer or a smart reference)
inline void MCParticle::addDecayMCVertex( MCVertex* value )                         {
  m_decayMCVertices.push_back(value); 
}
inline void MCParticle::addDecayMCVertex( SmartRef<MCVertex> value )                {
  m_decayMCVertices.push_back(value); 
}
*/

/// Serialize the object for writing
inline StreamBuffer& MCTrack::serialize( StreamBuffer& s ) const                 {
  ContainedObject::serialize(s);
// Ian  unsigned char u = (m_oscillationFlag) ? 1 : 0;
  return s
    << m_fourMomentum
    << m_particleID; // Ian Note need to get rid of the semicolon when adding junk
//    << u
//    << m_originMCVertex(this)
//    << m_decayMCVertices(this);
}


/// Serialize the object for reading
inline StreamBuffer& MCTrack::serialize( StreamBuffer& s )                       {
  ContainedObject::serialize(s);
// Ian  unsigned char u;
  s >> m_fourMomentum
    >> m_particleID;// Ian Note need to get rid of the semicolon when adding junk
//    >> u
//    >> m_originMCVertex(this)
//    >> m_decayMCVertices(this);
//  m_oscillationFlag = (u) ? true : false;
  return s;
}

/* 
/// Fill the ASCII output stream
inline std::ostream& MCParticle::fillStream( std::ostream& s ) const                {
  s << "class MCParticle :"
    << "\n    4-momentum (px, py, pz, E) = ( "
    << LHCbEventFloatFormat( LHCbEvent::width, LHCbEvent::precision )
    << m_fourMomentum.px() << ", "
    << LHCbEventFloatFormat( LHCbEvent::width, LHCbEvent::precision )
    << m_fourMomentum.py() << ", "
    << LHCbEventFloatFormat( LHCbEvent::width, LHCbEvent::precision )
    << m_fourMomentum.pz() << ", "
    << LHCbEventFloatFormat( LHCbEvent::width, LHCbEvent::precision )
    << m_fourMomentum.e()  << " )"
    << "\n    Particle ID                = " << m_particleID;
  if( m_oscillationFlag ) {
    s << m_oscillationFlag << "\n    Oscilation flag            = " << " true";
  }
  return s;
} */


// Definition of all container types of MCParticle
template <class TYPE> class ObjectVector;
typedef ObjectVector<MCTrack>     MCTrackVector;
template <class TYPE> class ObjectList;
typedef ObjectList<MCTrack>       MCTrackList;


#endif    // _H_GlastHitsEvt_MCTrack_v0_


