// $Header$
#ifndef LHCBEVENT_RAWEVENT_H
#define LHCBEVENT_RAWEVENT_H 1


// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "Gaudi/Kernel/DataObject.h"
#include "GlastEvent/Utilities/TriggerPattern.h"
#include "GlastEvent/TopLevel/Definitions.h"


// Externals 
extern const CLID& CLID_RawEvent;


//------------------------------------------------------------------------------
//
// ClassName:   RawEvent
//  
// Description: Essential information of the raw event
//              It can be identified by "/Event/Raw"
//
//              It contains:
//              - flag, if comming from Monte Carlo
//              - error status
//              - high voltage map
//              - trigger map
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------

/*!

Essential information of the raw event
It can be identified by "/Event/Raw"

It contains:
- flag, if comming from Monte Carlo
- error status
- high voltage map
- trigger map

 */

class RawEvent : public DataObject                                             {

public:
  /// Constructors
  RawEvent(const char* name = "RawEvent")
    : DataObject(name), 
      m_fromMC(false),
      m_errorStatus(0),
      m_highVoltageMask(0)                                                   { }
  /// Destructor
  virtual ~RawEvent()                                                        { }

  /// Retrieve reference to class definition structure
  virtual const CLID& clID() const               { return RawEvent::classID(); }
  static const CLID& classID()                         { return CLID_RawEvent; }

  ///  Retrieve flag of origine
  bool fromMC () const                                                         {
    return m_fromMC;
  }
  ///  Update flag of origine
  void setFromMC (bool value)                                                  {
    m_fromMC = value;
  }

  ///  Retrieve error status
  long errorStatus () const                                                    {
    return m_errorStatus;
  }
  ///  Update error status
  void setErrorStatus (long value)                                             {
    m_errorStatus = value;
  }

  ///  Retrieve high voltage status word
  long highVoltageMask () const                                                {
    return m_highVoltageMask;
  }
  ///  Update high voltage status word
  void setHighVoltageMask (long value)                                         {
    m_highVoltageMask = value;
  }

  ///  Retrieve trigger pattern
  const TriggerPattern& triggerPattern () const                                {
    return m_triggerPattern;
  }
  ///  Update trigger pattern
  void setTriggerPattern (const TriggerPattern& value)                         {
    m_triggerPattern = value;
  }

  /// Serialize the object for writing
  virtual StreamBuffer& serialize( StreamBuffer& s ) const;
  /// Serialize the object for reading
  virtual StreamBuffer& serialize( StreamBuffer& s );

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const RawEvent& obj )     {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  virtual std::ostream& fillStream( std::ostream& s ) const;

private: 
  /// Flag of origin
  bool              m_fromMC;
  /// Error status
  long              m_errorStatus;
  /// High voltage status word
  long              m_highVoltageMask;
  /// Trigger pattern
  TriggerPattern    m_triggerPattern;

};


//
// Inline code must be outside the class definition
//


/// Serialize the object for writing
inline StreamBuffer& RawEvent::serialize( StreamBuffer& s ) const              {
  DataObject::serialize(s);
  //unsigned char u = (m_fromMC) ? 1 : 0;
  return s;
//    << u
//    << m_errorStatus
//    << m_highVoltageMask
//    << m_triggerPattern;
}


/// Serialize the object for reading
inline StreamBuffer& RawEvent::serialize( StreamBuffer& s )                    {
  DataObject::serialize(s);
  unsigned char u;
 // s >> u
 //   >> m_errorStatus
 //   >> m_highVoltageMask
 //   >> m_triggerPattern;
 // m_fromMC = (u) ? true : false;
  return s;
}


/// Fill the output stream (ASCII)
inline std::ostream& RawEvent::fillStream( std::ostream& s ) const             {
  s << "class RawEvent :"
    << "\n    Flag of origin    = ";
  if( m_fromMC ) {
    s << " true";
  }
  else {
    s << "false";
  }
  s << "\n    Error status      = "
    << GlastEventField( GlastEvent::field12 )
    << m_errorStatus
    << "\n    High voltage mask = "
    << GlastEventField( GlastEvent::field12 )
    << m_highVoltageMask
    << "\n    Present triggers  = " << m_triggerPattern;
  return s;
}


#endif  // LHCBEVENT_RAWEVENT_H
