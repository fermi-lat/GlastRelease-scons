// $Header$
#ifndef LHCBEVENT_TIMESTAMP_H
#define LHCBEVENT_TIMESTAMP_H 1


// Include files
#include <iostream>
#include "Gaudi/Kernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/Definitions.h"


//------------------------------------------------------------------------------
//
// ClassName:   TimeStamp
//  
// Description: Placeholder for the time stamp
//              (decision about precision and range has to be made)
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class TimeStamp                                                                {

public:

  /// Constructors
  TimeStamp()
    : m_time(0)                                                              { }
  TimeStamp( long t )
    : m_time(t)                                                              { }
  /// Destructor
  ~TimeStamp()                                                               { }

  /// Retrieve time
  long time() const                                                            {
    return m_time;
  }
  /// Update time 
  void setTime( long value )                                                   {
    m_time = value;
  }

  /// Serialize the object for writing
  friend StreamBuffer& operator<< ( StreamBuffer& s, const TimeStamp& obj )    {
    return s << obj.m_time;
  }
  /// Serialize the object for reading
  friend StreamBuffer& operator>> ( StreamBuffer& s, TimeStamp& obj )          {
    return s >> obj.m_time;
  }

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const TimeStamp& obj )    {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  std::ostream& fillStream( std::ostream& s ) const                            {
    return s << "class TimeStamp : "
      << GlastEventField( GlastEvent::field12 )
      << m_time;
  }

private:

  /// Time
  long m_time;

};


#endif    // LHCBEVENT_TIMESTAMP_H
