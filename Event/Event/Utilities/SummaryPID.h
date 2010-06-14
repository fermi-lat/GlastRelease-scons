// $Header$
#ifndef LHCBEVENT_SUMMARYPID_H
#define LHCBEVENT_SUMMARYPID_H 1


// Include files
#include <iostream>
#include "Gaudi/Kernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/Definitions.h"


//------------------------------------------------------------------------------
//
// ClassName:   SummaryPId
// 
// Description: Class containg summary of all possible compatible
//              particleID with a track
//              For the moment it is a placeholder containg 2 
//              integer values, where coded words for possible particle
//              and antiparticleID are stored
//
// Author:      G.Corti
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class SummaryPID                                                               {

public:

  // Constructors
  SummaryPID( long pValue, long antipValue )
    : m_pIDCode(pValue),
      m_antiPIDCode(antipValue)                                              { }
  SummaryPID()
    : m_pIDCode(0),
      m_antiPIDCode(0)                                                       { }
  // Destructor
  ~SummaryPID()                                                              { }

  // Retrieve first element of summary particle ID
  // at the moment coded word containing all comaptible particle IDs
  long pIDCode() const                                                       {
    return m_pIDCode;
  }

  /// Update first element of summary particle ID
  void setPIDCode( long value )                                              {
    m_pIDCode = value;
  }

  // Retrieve second element of summary particle Id
  // at the moment coded word containing all comaptible particle IDs
  long antiPIDCode() const                                                   {
    return m_antiPIDCode;
  }

  /// Update second element of summary particle ID
  void setAntiPIDCode( long value )                                          {
    m_antiPIDCode = value;
  }

  /// Serialize the object for writing
  friend StreamBuffer& operator<< ( StreamBuffer& s, const SummaryPID& obj ) {
    return s
      << obj.m_pIDCode
      << obj.m_antiPIDCode;
  }

  /// Serialize the object for reading
  friend StreamBuffer& operator>> ( StreamBuffer& s, SummaryPID& obj )       {
    return s
      >> obj.m_pIDCode
      >> obj.m_antiPIDCode;
  }

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const SummaryPID& obj ) {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  std::ostream& fillStream( std::ostream& s ) const                          {
    return s
      << "class SummaryPID : ( "
      << GlastEventField( GlastEvent::field5 )
      << m_pIDCode << ", "
      << GlastEventField( GlastEvent::field5 )
      << m_antiPIDCode << " )";
  }

private:

  long m_pIDCode;                 // Particle ID Summary number
  long m_antiPIDCode;             // Anti Particle ID Summary number

};


#endif    // LHCBEVENT_PARTICLEID_H
