// $Header$
#ifndef LHCBEVENT_PARTICLEID_H
#define LHCBEVENT_PARTICLEID_H 1


// Include files
#include <iostream>
#include "Gaudi/Kernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/Definitions.h"


//------------------------------------------------------------------------------
//
// ClassName:   ParticleID
//  
// Description: Particle identifier
//              (corresponding to the particle identification table)
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class ParticleID                                                               {

public:

  /// Constructors
  ParticleID( long id )
    : m_id(id)                                                               { }
  ParticleID()
    : m_id(0)                                                                { }
  /// Destructor
  ~ParticleID()                                                              { }

  /// Retrieve particle identifiaction
  long id () const                                                             {
    return m_id;
  }
  /// Update particle identifiaction
  void setID (long value)                                                      {
    m_id = value;
  }

  /// Serialize the object for writing
  friend StreamBuffer& operator<< ( StreamBuffer& s, const ParticleID& obj )   {
    return s << obj.m_id;
  }
  /// Serialize the object for reading
  friend StreamBuffer& operator>> ( StreamBuffer& s, ParticleID& obj )         {
    return s >> obj.m_id;
  }

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const ParticleID& obj )   {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  std::ostream& fillStream( std::ostream& s ) const                            {
    return s << "class ParticleID : "
      << GlastEventField( GlastEvent::field4 )
      << m_id;
  }

private:

  /// Particle identifiaction number
  long m_id;

};


#endif    // LHCBEVENT_PARTICLEID_H
