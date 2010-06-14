// $Header$
#ifndef LHCBEVENT_ANALEVENT_H
#define LHCBEVENT_ANALEVENT_H 1


// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/DataObject.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "GlastEvent/Utilities/Classification.h"
#include "GlastEvent/TopLevel/Definitions.h"


// Externals 
extern const CLID& CLID_AnalEvent;


//------------------------------------------------------------------------------
//
// ClassName:   AnalEvent
//  
// Description: Essential information of the analysis event
//              It can be identified by "/Event/Anal"
//
//              It contains:
//              - event classification
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class AnalEvent : public DataObject                                            {

public:
  /// Constructors
  AnalEvent(const char* name = "AnalEvent")
    : DataObject(name)                                                       { }
  /// Destructor
  virtual ~AnalEvent()                                                       { }

  /// Retrieve reference to class definition structure
  virtual const CLID& clID() const              { return AnalEvent::classID(); }
  static const CLID& classID()                        { return CLID_AnalEvent; }

  /// Retrieve classificiation
  const Classification& classification() const                                 {
    return m_classification;
  }
  /// Update classification
  void setClassification (const Classification& value)                         {
    m_classification = value;
  }

  /// Serialize the object for writing
  virtual StreamBuffer& serialize( StreamBuffer& s ) const;
  /// Serialize the object for reading
  virtual StreamBuffer& serialize( StreamBuffer& s );

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const AnalEvent& obj )    {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  virtual std::ostream& fillStream( std::ostream& s ) const;

private:
  /// Event classification
  Classification    m_classification;

};


//
// Inline code must be outside the class definition
//


/// Serialize the object for writing
inline StreamBuffer& AnalEvent::serialize( StreamBuffer& s ) const             {
  DataObject::serialize(s);
  return s << m_classification;
}


/// Serialize the object for reading
inline StreamBuffer& AnalEvent::serialize( StreamBuffer& s )                   {
  DataObject::serialize(s);
  return s >> m_classification;
}


/// Fill the output stream (ASCII)
inline std::ostream& AnalEvent::fillStream( std::ostream& s ) const            {
  return s
    << "class AnalEvent :"
    << "\n    Event classification = "
    << GlastEventField( GlastEvent::field12 )
    << m_classification;
}


#endif      // LHCBEVENT_ANALEVENT_H
