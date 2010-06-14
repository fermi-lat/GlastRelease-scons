// $Header$
#ifndef LHCBEVENT_CLASSIFICATION_H
#define LHCBEVENT_CLASSIFICATION_H 1


// Include files
#include <iostream>
#include "Gaudi/Kernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/Definitions.h"


//------------------------------------------------------------------------------
//
// ClassName:   Classification
//  
// Description: Placeholder for (sub)event classification
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class Classification                                                           {

public:

  /// Constructors
  Classification (long eventType)
    : m_eventType(eventType)                                                 { }
  Classification()
    : m_eventType(0)                                                         { }
  /// Destructor
  ~Classification()                                                          { }

  /// Retrieve event type information
  long eventType() const                                                       {
    return m_eventType;
  }
  /// Update event type information
  void setEventType( long value )                                              {
    m_eventType = value;
  }

  /// Serialize the object for writing
  friend StreamBuffer& operator<< ( StreamBuffer& s,
                                    const Classification& obj )                {
    return s << obj.m_eventType;
  }
  /// Serialize the object for reading
  friend StreamBuffer& operator>> ( StreamBuffer& s,
                                    Classification& obj )                      {
    return s >> obj.m_eventType;
  }
  
  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s,
                                    const Classification& obj )                {
    return obj.fillStream(s);
  }  
  /// Fill the output stream (ASCII)
  std::ostream& fillStream( std::ostream& s ) const                            {
    return s << "class Classification : "
      << GlastEventField( GlastEvent::field4 )
      << m_eventType;
  }

private:

  /// Event type info
  long m_eventType;

};


#endif    // LHCBEVENT_CLASSIFICATION_H
