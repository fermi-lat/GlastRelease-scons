// $Header$
#ifndef LHCBEVENT_RECEVENT_H
#define LHCBEVENT_RECEVENT_H 1


// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "Gaudi/Kernel/DataObject.h"
#include "GlastEvent/TopLevel/Definitions.h"


// Externals 
extern const CLID& CLID_RecEvent;


//------------------------------------------------------------------------------
//
// ClassName:   RecEvent
//  
// Description: Essential information of the reconstructed event
//              It can be identified by "/Event/Rec"
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class RecEvent : public DataObject                                             {

public:
  /// Constructors
  RecEvent( const char* name = "RecEvent" )
    : DataObject(name)                                                       { }
  /// Destructor
  virtual ~RecEvent()                                                        { }

  /// Retrieve reference to class definition structure
  virtual const CLID& clID() const               { return RecEvent::classID(); }
  static const CLID& classID()                         { return CLID_RecEvent; }

  /// Serialize the object for writing
  virtual StreamBuffer& serialize( StreamBuffer& s ) const;
  /// Serialize the object for reading
  virtual StreamBuffer& serialize( StreamBuffer& s );

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const RecEvent& obj )     {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  virtual std::ostream& fillStream( std::ostream& s ) const;

private:
  /// Dummy data member
  double                        m_dummy;

};


//
// Inline code must be outside the class definition
//


/// Serialize the object for writing
inline StreamBuffer& RecEvent::serialize( StreamBuffer& s ) const              {
  DataObject::serialize(s);
  return s << m_dummy;
}


/// Serialize the object for reading
inline StreamBuffer& RecEvent::serialize( StreamBuffer& s )                    {
  DataObject::serialize(s);
  return s >> m_dummy;
}


/// Fill the output stream (ASCII)
inline std::ostream& RecEvent::fillStream( std::ostream& s ) const             {
  return s
    << "class RecEvent :"
    << "\n    Dummy     = "
    << GlastEventField( GlastEvent::field12 )
    << m_dummy;
}


#endif    // LHCBEVENT_RECEVENT_H
