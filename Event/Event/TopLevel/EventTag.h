// $Header$
#ifndef LHCBEVENT_EVENTTAG_H
#define LHCBEVENT_EVENTTAG_H 1


// Include files
#include <iostream>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/SmartRef.h"
#include "Event/Utilities/TriggerPattern.h"
#include "Event/Utilities/Classification.h"
#include "Event/TopLevel/Definitions.h"


// Forward declarations
class Event;


// Externals 
extern const CLID& CLID_EventTag;


//------------------------------------------------------------------------------
//
// ClassName:   EventTag
//  
// Description: Essential event tag information
//
//              It contains:
//              - trigger pattern
//              - event classification
//              - reference to the mother event
//
// Author:      Pavel Binko
// Changes:     M.Frank 04/10/1999 : Proper use of SmartRefs and SmartRefVectors
//              P.Binko 19/10/1999 : Proper accessors of smart references,
//                                   Formating of ASCII output
//
//------------------------------------------------------------------------------

//! Essential event tag information
/*!
It contains:
- trigger pattern
- event classification
- reference to the mother event

*/



class EventTag : public DataObject                                             {

public:
  /// Constructor
  EventTag(Event* e, TriggerPattern p, Classification c)
    : m_event(e),
      m_triggerPattern(p),
      m_classification(c)                                                    { }
  /// Destructor
  ~EventTag()                                                                { }

  /// Retrieve reference to class definition structure
  virtual const CLID& clID() const               { return EventTag::classID(); }
  static const CLID& classID()                         { return CLID_EventTag; }

  /// Retrieve event trigger pattern
  const TriggerPattern& triggerPattern () const                                {
    return m_triggerPattern;
  }
  /// Update event trigger pattern
  void setTriggerPattern (const TriggerPattern& value)                         {
    m_triggerPattern = value;
  }
  
  /// Retrieve event classification
  const Classification& classification () const                                {
    return m_classification;
  }
  /// Update event classification
  void setClassification (Classification value)                                {
    m_classification = value;
  }
  
  /// Retrieve pointer to event structure (const or non-const)
  const Event* EventTag::event() const;
        Event* EventTag::event();
  /// Update pointer to event structure (by a C++ pointer or a smart reference)
  void EventTag::setEvent( Event* value );
  void EventTag::setEvent( SmartRef<Event> value );

  /// Serialize the object for writing
  virtual StreamBuffer& serialize( StreamBuffer& s ) const;
  /// Serialize the object for reading
  virtual StreamBuffer& serialize( StreamBuffer& s );

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const EventTag& obj )     {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  std::ostream& fillStream( std::ostream& s ) const;

private:
  /// Trigger pattern
  TriggerPattern   m_triggerPattern;
  /// Event classification
  Classification   m_classification;
  /// Reference to the event
  SmartRef<Event>  m_event;

};


//
// Inline code must be outside the class definition
//
#include "GlastEvent/TopLevel/Event.h"


/// Retrieve pointer to event structure (const or non-const)
inline const Event* EventTag::event() const                                    {
  return m_event; 
}
inline       Event* EventTag::event()                                          {
  return m_event; 
}
/// Update pointer to event structure (by a C++ pointer or a smart reference)
inline void EventTag::setEvent( Event* value )                                 {
  m_event = value; 
}
inline void EventTag::setEvent( SmartRef<Event> value )                        {
  m_event = value; 
}


/// Serialize the object for writing
inline StreamBuffer& EventTag::serialize( StreamBuffer& s ) const              {
  DataObject::serialize(s);
  return s
    << m_triggerPattern
    << m_classification
    << m_event(this);
}


/// Serialize the object for reading
inline StreamBuffer& EventTag::serialize( StreamBuffer& s )                    {
  DataObject::serialize(s);
  return s
    >> m_triggerPattern
    >> m_classification
    >> m_event(this);
}


/// Fill the output stream (ASCII)
inline std::ostream& EventTag::fillStream( std::ostream& s ) const             {
  return s
    << "class EventTag :"
    << "\n    Present triggers     = " << m_triggerPattern
    << "\n    Event classification = " << m_classification;
}


#endif      // LHCBEVENT_EVENTTAG_H
