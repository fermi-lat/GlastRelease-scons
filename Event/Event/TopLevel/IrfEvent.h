#ifndef LHCBEVENT_IrfEvent_H
#define LHCBEVENT_IrfEvent_H 1


// Include files
#include <iostream>
#include <vector>
#include <algorithm>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/DataObject.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/Definitions.h"


// Externals
extern const CLID& CLID_IrfEvent;


//------------------------------------------------------------------------------
//
// ClassName:   IrfEvent
//  
// Description: Essential information of the IRF event
//              It can be identified by "/Event/Irf"
//
//              It contains:
//
//
//------------------------------------------------------------------------------

/*!

Essential information of the IRF event
It can be identified by "/Event/Irf"

*/

class IrfEvent : public DataObject                                              {

public:
  /// Constructors
    // HMA removed 2nd parameter pileUp - not sure what it was doing...
  IrfEvent( const char* name = "IrfEvent" )
      : DataObject(name) {}
  /// Destructor
  virtual ~IrfEvent()                                                         { }

  /// Retrieve reference to class definition structure
  virtual const CLID& clID() const                { return IrfEvent::classID(); }
  static const CLID& classID()                          { return CLID_IrfEvent; }

  /// Clone operator
  IrfEvent& operator=(const IrfEvent& copy)                                      {
    return *this;
  }

  /// Serialize the object for writing
  virtual StreamBuffer& serialize( StreamBuffer& s ) const;
  /// Serialize the object for reading
  virtual StreamBuffer& serialize( StreamBuffer& s );

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const IrfEvent& obj )      {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  virtual std::ostream& fillStream( std::ostream& s ) const;

private:

};


//
// Inline code must be outside the class definition
//

/// Serialize the object for writing
inline StreamBuffer& IrfEvent::serialize( StreamBuffer& s ) const               {
  DataObject::serialize(s);
  // if there was something to serialize we'd do it here
  return s;
}


/// Serialize the object for reading
inline StreamBuffer& IrfEvent::serialize( StreamBuffer& s )                     {
  DataObject::serialize(s);
// if there were something to serialize we'd do it here
  return s;
}


/// Fill the ASCII output stream
inline std::ostream& IrfEvent::fillStream( std::ostream& s ) const              {
  s << "class IrfEvent :\n";

  return s;
}


#endif    // GLASTEVENT_IrfEvent_H
