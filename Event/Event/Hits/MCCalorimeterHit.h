// $Header$
#ifndef MCCalorimeterHit_H
#define MCCalorimeterHit_H 1


// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "Gaudi/Kernel/SmartRefVector.h"
#include "GlastEvent/TopLevel/Definitions.h"
#include "GlastEvent/Utilities/CellID.h"
#include "GlastEvent/Hits/MCTrack.h"
// Include all LHCb container types here
//   to simplify inlude statements in algorithms
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"


// Forward declarations
//class MCParticle;


// Externals 
extern const CLID& CLID_MCCalorimeterHit;


//------------------------------------------------------------------------------
//
// ClassName:   MCCalorimeterHit
//  
// Description: Monte Carlo calorimeter hit
//              Definition for hits from HCal hits  (SICB HCEL and HCMT banks)
//                                       Ecal hits  (SICB ECEL and ECMT banks)
//                                  preshower hits  (SICB ECPC and HCMT banks)
//
//              Can evolve into a base class of ECal, HCal and preshower hits
//
// Author:      Pavel Binko
//
// Changes:     M.Frank 04/10/1999 : Proper use of SmartRefs and SmartRefVectors
//              P.Binko 19/10/1999 : Proper accessors of smart references,
//                                   Formating of ASCII output
//
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//
// GlastName:	MCCalorimeterHits
//
// Changes:		Modified to use the MCTrack class rather than the MCParticle 
//				class from LHCbEvent
//
//
//------------------------------------------------------------------------------

class MCCalorimeterHit : virtual public ContainedObject {

public:
  /// Constructors
  MCCalorimeterHit()                                    { }
  /// Destructor
  virtual ~MCCalorimeterHit()                           { }

  /// Retrieve pointer to class defininition structure
  virtual const CLID& clID() const   { return MCCalorimeterHit::classID(); }
  static const CLID& classID()       { return CLID_MCCalorimeterHit; }

  /// Retrieve energy
  double energy() const                                 { return m_energy; }
  /// Update energy
  void setEnergy( double value )                        { m_energy = value; }
  /// Retrieve cell identifier
  const CellID cellID() const                           { return m_cellID; }
  /// Update cell identifier
  void setCellID( CellID value )                        { m_cellID = value; }

  /// Retrieve pointer to vector of MCParticles (const or non-const)
  const SmartRefVector<MCTrack>& mcTracks() const;
        SmartRefVector<MCTrack>& mcTracks();
  /// Update all MCParticles
  void setMCTracks( const SmartRefVector<MCTrack>& value );
  /// Remove all MCParticles
  void removeMCTracks();
  /// Add single MCParticle to SmartRefVector<MCParticle>
  ///   (by a C++ pointer or a smart reference)
  void addMCTrack( MCTrack* value );
  void addMCTrack( SmartRef<MCTrack> value );

  /// Serialize the object for reading
  virtual StreamBuffer& serialize( StreamBuffer& s );
  /// Serialize the object for writing
  virtual StreamBuffer& serialize( StreamBuffer& s ) const;
  /// Fill the ASCII output stream
  virtual std::ostream& fillStream( std::ostream& s ) const;

  
  // response from the minus end of the log
  float minusEnd () const             { return m_minusEnd; }
  void setMinusEnd(float value)       { m_minusEnd = value; }

  // reponse from the plus end of the log
  float plusEnd() const               { return m_plusEnd; }
  void setPlusEnd(float value)        {m_plusEnd = value;}

  
  unsigned int layer () const         { return m_layer; }
  void setLayer(unsigned int value) { m_layer = value; }

  /// Log index number
  unsigned int index () const           { return m_index; }
  void setId (long value)            { m_id = value; }



private:
  /// Energy
  double                        m_energy;
  /// Cell identifier
  CellID                        m_cellID;
  /// Vector of MCParticles
  SmartRefVector<MCTrack>    m_mcTracks;
  
  float                m_minusEnd; // response from the diode on the "Minus" end of the log
  float                m_plusEnd;  // response from the diode on the "Plus" end of the log
  
  unsigned int         m_index;    // index of this crystal within a layer
  unsigned int         m_layer;    // layer number
  unsigned int         m_id;	   // Log index number


};


//
// Inline code must be outside the class definition
//
#include "GlastEvent/Hits/MCTrack.h"


/// Retrieve pointer to vector of MCParticles (const or non-const)
inline const SmartRefVector<MCTrack>& MCCalorimeterHit::mcTracks() const            {
  return m_mcTracks;
}
inline       SmartRefVector<MCTrack>& MCCalorimeterHit::mcTracks()                  {
  return m_mcTracks;
}
/// Update all MCParticles
inline void MCCalorimeterHit::setMCTracks( const SmartRefVector<MCTrack>& value )   {
  m_mcTracks = value;
}
/// Remove all MCParticles
inline void MCCalorimeterHit::removeMCTracks()                                         {
  m_mcTracks.clear();
}
/// Add single MCParticle to SmartRefVector<MCParticle>
///   (by a C++ pointer or a smart reference)
inline void MCCalorimeterHit::addMCTrack( MCTrack* value )                          {
  m_mcTracks.push_back(value);
}
inline void MCCalorimeterHit::addMCTrack( SmartRef<MCTrack> value )                 {
  m_mcTracks.push_back(value);
}


/// Serialize the object for writing
inline StreamBuffer& MCCalorimeterHit::serialize( StreamBuffer& s ) const                 {
  ContainedObject::serialize(s);
  return s
    << m_cellID
    << m_energy
    << m_mcTracks(this);
}


/// Serialize the object for reading
inline StreamBuffer& MCCalorimeterHit::serialize( StreamBuffer& s )                       {
  ContainedObject::serialize(s);
  return s
    >> m_cellID
    >> m_energy
    >> m_mcTracks(this);
}


/// Fill the ASCII output stream
inline std::ostream& MCCalorimeterHit::fillStream( std::ostream& s ) const                {
  return s
    << "class MCCalorimeterHit :"
    << "\n    Energy    = "
//    << LHCbEventFloatFormat( LHCbEvent::width, LHCbEvent::precision )
    << m_energy
    << "\n    Cell ID   = " << m_cellID;
}

  
// Definition of all container types of MCCalorimeterHit
template <class TYPE> class ObjectVector;
typedef ObjectVector<MCCalorimeterHit>     MCCalorimeterHitVector;
template <class TYPE> class ObjectList;
typedef ObjectList<MCCalorimeterHit>       MCCalorimeterHitList;


#endif // GLASTEVENT_MCCALORIMETERHIT_H
