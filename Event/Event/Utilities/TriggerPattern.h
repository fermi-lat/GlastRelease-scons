// $Header$
#ifndef LHCBEVENT_TRIGGERPATTERN_H
#define LHCBEVENT_TRIGGERPATTERN_H 1


// Include files
#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"


//------------------------------------------------------------------------------
//
// ClassName:   TriggerPattern
//  
// Description: Placeholder for trigger pattern
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class TriggerPattern                                                           {

public:

  /// Constructors
  TriggerPattern()
    : m_pattern(0)                                                           { }
  TriggerPattern( long pattern )
    : m_pattern(pattern)                                                     { }
  /// Destructor
  ~TriggerPattern()                                                          { }

  /// Retrieve trigger bits
  long pattern() const                                                         {
    return m_pattern;
  }
  /// Update trigger bits
  void setPattern( long value )                                                {
    m_pattern = value;
  }
  /// Set single bit
  void set( int bit )                                                          {
    if ( bit>=0 && bit<32 ) {
      m_pattern |= (1<<bit);
    }
  }
  /// Check specified bit
  bool isSet( int bit ) const                                                  {
    return (bit>=0 && bit<32) ? (1==((m_pattern&(1<<bit)) >> bit)) : false;
  }

  /// Serialize the object for writing
  friend StreamBuffer& operator<< ( StreamBuffer& s,
                                    const TriggerPattern& obj )                {
    return s << obj.m_pattern;
  }
  /// Serialize the object for reading
  friend StreamBuffer& operator>> ( StreamBuffer& s,
                                    TriggerPattern& obj )                      {
    return s >> obj.m_pattern;
  }

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s,
                                    const TriggerPattern& obj )                {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  std::ostream& fillStream( std::ostream& s ) const                            {
    return s << "class TriggerPattern : "
      << GlastEventField( EventFormat::field12 )
      << m_pattern;
  }

private:

  /// Data member to hold trigger bits
  long m_pattern;

};


#endif  // LHCBEVENT_TRIGGERPATTERN_H
