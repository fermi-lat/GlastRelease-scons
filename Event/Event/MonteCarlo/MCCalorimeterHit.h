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
#include "GlastEvent/MonteCarlo/MCTrack.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"

/*! Represents the data from the CsI logs. Originally adapted from H. Arrighi
    MCACDHit class. The IRF parser used with this class is CsIDetector
*/


// Forward declarations
//class MCParticle;


//! Externals 
extern const CLID& CLID_MCCalorimeterHit;




//------------------------------------------------------------------------------
//
// GlastName:   MCCalorimeterHits
//
// Notes:       This class was reworked from the class of the same name in  
//              LHCbEvent
//
//
//------------------------------------------------------------------------------

class MCCalorimeterHit : virtual public ContainedObject {

public:
  //! Constructors
  MCCalorimeterHit()                                    { }
  //! Destructor
  virtual ~MCCalorimeterHit()                           { }

  //! Retrieve pointer to class defininition structure
  virtual const CLID& clID() const   { return MCCalorimeterHit::classID(); }
  static const CLID& classID()       { return CLID_MCCalorimeterHit; }

  //! Retrieve energy
  double energy() const                                 { return m_energy; }
  //! Update energy
  void setEnergy( double value )                        { m_energy = value; }
  //! Retrieve cell identifier
  const CellID cellID() const                           { return m_cellID; }
  //! Update cell identifier
  void setCellID( CellID value )                        { m_cellID = value; }

  //! Retrieve pointer to vector of MCTrack (const or non-const)
  const SmartRefVector<MCTrack>& mcTracks() const;
        SmartRefVector<MCTrack>& mcTracks();
  //! Update all MCParticles
  void setMCTracks( const SmartRefVector<MCTrack>& value );
  //! Remove all MCParticles
  void removeMCTracks();
  /*! Add single MCParticle to SmartRefVector<MCTrack>
     (by a C++ pointer or a smart reference)
  */
  void addMCTrack( MCTrack* value );
  void addMCTrack( SmartRef<MCTrack> value );

  //! Serialize the object for reading
  virtual StreamBuffer& serialize( StreamBuffer& s );
  //! Serialize the object for writing
  virtual StreamBuffer& serialize( StreamBuffer& s ) const;
  //! Fill the ASCII output stream
  virtual std::ostream& fillStream( std::ostream& s ) const;

  
  //! response from the right end of the log
  double rightResponse() const               { return m_rightEnd; }
  void  setRightResponse(float value)       { m_rightEnd = value; }

  //! reponse from the left end of the log
  double leftResponse() const                { return m_leftEnd; }
  void  setLeftResponse(float value)        {m_leftEnd = value;}

  //!  The Layer
  unsigned int layer () const               { return m_layer; }
  void setLayer(unsigned int value)         { m_layer = value; }

  //! Log index number
  unsigned int index () const               { return m_index; }
  void setId (long value)                   { m_id = value; }



private:
  //! Energy
  double                        m_energy;  
  double                        m_leftEnd;  // response from the diode on the "left" end of the log
  double                        m_rightEnd; // response from the diode on the "right" end of the log
  //! Cell identifier
  CellID                        m_cellID;
  //! Vector of MCParticles
  SmartRefVector<MCTrack>    m_mcTracks;
  

  
  unsigned int         m_index;    // index of this crystal within a layer
  unsigned int         m_layer;    // layer number
  unsigned int         m_id;	   // Log index number


};


//
// Inline code must be outside the class definition
//
#include "GlastEvent/MonteCarlo/MCTrack.h"


//! Retrieve pointer to vector of MCParticles (const or non-const)
inline const SmartRefVector<MCTrack>& MCCalorimeterHit::mcTracks() const            {
  return m_mcTracks;
}
inline       SmartRefVector<MCTrack>& MCCalorimeterHit::mcTracks()                  {
  return m_mcTracks;
}
/// Update all MCTracks
inline void MCCalorimeterHit::setMCTracks( const SmartRefVector<MCTrack>& value )   {
  m_mcTracks = value;
}
//! Remove all MCTracks
inline void MCCalorimeterHit::removeMCTracks()                                         {
  m_mcTracks.clear();
}
/*! Add single MCParticle to SmartRefVector<MCParticle>
   (by a C++ pointer or a smart reference)*/
inline void MCCalorimeterHit::addMCTrack( MCTrack* value )                          {
  m_mcTracks.push_back(value);
}
inline void MCCalorimeterHit::addMCTrack( SmartRef<MCTrack> value )                 {
  m_mcTracks.push_back(value);
}


//! Serialize the object for writing
inline StreamBuffer& MCCalorimeterHit::serialize( StreamBuffer& s ) const                 {
  ContainedObject::serialize(s);
  return s
//    << m_cellID
    << m_energy
    << m_leftEnd
    << m_rightEnd;
//    << m_mcTracks(this);
}


//! Serialize the object for reading
inline StreamBuffer& MCCalorimeterHit::serialize( StreamBuffer& s )                       {
  ContainedObject::serialize(s);
  return s
//    >> m_cellID
    >> m_energy
    >> m_leftEnd
    >> m_rightEnd;
//    >> m_mcTracks(this);
}


//! Fill the ASCII output stream
inline std::ostream& MCCalorimeterHit::fillStream( std::ostream& s ) const                {
  return s
    << "class MCCalorimeterHit :"
    << "\n    Energy    = "
    << m_energy
    << m_leftEnd
    << m_rightEnd
    << "\n    Cell ID   = " << m_cellID;
}

  
//! Definition of all container types of MCCalorimeterHit
template <class TYPE> class ObjectVector;
typedef ObjectVector<MCCalorimeterHit>     MCCalorimeterHitVector;
template <class TYPE> class ObjectList;
typedef ObjectList<MCCalorimeterHit>       MCCalorimeterHitList;


#endif // GLASTEVENT_MCCALORIMETERHIT_H
