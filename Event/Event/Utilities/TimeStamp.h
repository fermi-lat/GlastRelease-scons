#ifndef EVENT_TIMESTAMP_H
#define EVENT_TIMESTAMP_H 1

#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"

/** @class TimeStamp
* @brief encapsulate the time.
*
* Elapsed? absolute? Currently a double, in units of seconds.
*
* $Header$
*/

class TimeStamp                                                                {

public:

  TimeStamp()
    : m_time(0)                                                              { }
  TimeStamp( double t )
    : m_time(t)                                                              { }

  ~TimeStamp()                                                               { }

  /// Retrieve time
  double time() const                                                            {
    return m_time;
  }
  /// Update time 
  void setTime( double value )                                                   {
    m_time = value;
  }

  operator double()const { return time(); }

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
      << EventField( EventFormat::field12 )
      << m_time;
  }

private:

  /// Time
  double m_time;

};


#endif    // LHCBEVENT_TIMESTAMP_H
