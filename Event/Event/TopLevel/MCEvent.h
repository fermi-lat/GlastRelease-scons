// $Header$
#ifndef LHCBEVENT_MCEVENT_H
#define LHCBEVENT_MCEVENT_H 1


// Include files
#include <iostream>
#include <vector>
#include <algorithm>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/DataObject.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/SubMCEvent.h"
#include "GlastEvent/TopLevel/Definitions.h"


// Externals
extern const CLID& CLID_MCEvent;


//------------------------------------------------------------------------------
//
// ClassName:   MCEvent
//  
// Description: Essential information of the Monte Carlo event
//              It can be identified by "/Event/MC"
//
//              It contains:
//              - pile up
//              - vector of Monte Carlo sub-events
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------

/*!

Essential information of the Monte Carlo event
It can be identified by "/Event/MC"

It contains:
- pile up
- vector of Monte Carlo sub-events
*/

class MCEvent : public DataObject                                              {

public:
  /// Constructors
    // HMA removed 2nd parameter pileUp - not sure what it was doing...
  MCEvent( const char* name = "MCEvent" )
      : DataObject(name) {}
    //m_pileUp(pileUp)                                                         { }
  /// Destructor
  virtual ~MCEvent()                                                         { }

  /// Retrieve reference to class definition structure
  virtual const CLID& clID() const                { return MCEvent::classID(); }
  static const CLID& classID()                          { return CLID_MCEvent; }

  /// Clone operator
  MCEvent& operator=(const MCEvent& copy)                                      {
    //m_pileUp       = copy.pileUp();
    m_subMCEvents  = copy.subMCEvent();
    return *this;
  }

  /*
  /// Retrieve pileUp
  long pileUp () const                                                         {
    return m_pileUp;
  }
  /// Update pileUp
  void setPileUp (long value)                                                  {
    m_pileUp = value;
  }
*/
  /// Retrieve pointer to vector of MC sub event entries (const or non-const)
  const std::vector<SubMCEvent*>& subMCEvent() const;
        std::vector<SubMCEvent*>& subMCEvent();
  /// Update all MC sub event entries
  void setSubMCEvents( const std::vector<SubMCEvent*>& value );
  /// Remove all MC sub event entries
  void removeSubMCEvents();
  /// Add single MC sub event entry to vector
  void addSubMCEvent( SubMCEvent* value );
  /// Remove single MC sub event entry from vector
  ///   This function does not delete the object "value" point to
  void removeSubMCEvent( SubMCEvent* value );

  /// Serialize the object for writing
  virtual StreamBuffer& serialize( StreamBuffer& s ) const;
  /// Serialize the object for reading
  virtual StreamBuffer& serialize( StreamBuffer& s );

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const MCEvent& obj )      {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  virtual std::ostream& fillStream( std::ostream& s ) const;

private:
  /// Pile up flags
  long                      m_pileUp;
  /// Vector of Monte Carlo sub events
  std::vector<SubMCEvent*>  m_subMCEvents;

};


//
// Inline code must be outside the class definition
//
//#include "GlastEvent/MonteCarlo/MCVertex.h"  HMA commented this out until we have MCVertex.h


/// Retrieve pointer to vector of MC sub event entries (const or non-const)
inline const std::vector<SubMCEvent*>& MCEvent::subMCEvent() const             { 
  return m_subMCEvents; 
}
inline       std::vector<SubMCEvent*>& MCEvent::subMCEvent()                   { 
  return m_subMCEvents; 
}
/// Update all MC sub event entries
inline void MCEvent::setSubMCEvents( const std::vector<SubMCEvent*>& value )   { 
  m_subMCEvents = value; 
}
/// Remove all MC sub event entries
inline void MCEvent::removeSubMCEvents()                                       { 
  m_subMCEvents.clear(); 
}
/// Add single MC sub event entry to vector
inline void MCEvent::addSubMCEvent( SubMCEvent* value )                        {
  m_subMCEvents.push_back(value); 
}
/// Remove single MC sub event entry from vector
///   This function does not delete the object "value" point to
inline void MCEvent::removeSubMCEvent( SubMCEvent* value )                     {
  std::vector<SubMCEvent*>::iterator i = 
    std::remove( m_subMCEvents.begin(), m_subMCEvents.end(), value );
  // If something removed, erase what has left at the of our vector
  if ( i != m_subMCEvents.end() ) {
    m_subMCEvents.erase( i, m_subMCEvents.end() );
  }
}


/// Serialize the object for writing
inline StreamBuffer& MCEvent::serialize( StreamBuffer& s ) const               {
  DataObject::serialize(s);
  s //<< m_pileUp
    << m_subMCEvents.size();
  std::vector<SubMCEvent*>::const_iterator iter;
  for( iter = m_subMCEvents.begin(); iter != m_subMCEvents.end(); iter++ ) {
    s << **iter;
  }
  return s;
}


/// Serialize the object for reading
inline StreamBuffer& MCEvent::serialize( StreamBuffer& s )                     {
  DataObject::serialize(s);
  std::vector<SubMCEvent*>::size_type    siz;
  SubMCEvent*                            pSubMCEvent;
  s //>> m_pileUp
    >> siz;
  for( long i = 0; i < (long)siz; i++ ) {
    pSubMCEvent = new SubMCEvent;
    s >> *pSubMCEvent;
    m_subMCEvents.push_back(pSubMCEvent);
  }
  return s;
}


/// Fill the ASCII output stream
inline std::ostream& MCEvent::fillStream( std::ostream& s ) const              {
  s << "class MCEvent :\n"
    << "    Pile-up = "
    << GlastEventField( GlastEvent::field12 );
    //<< m_pileUp;
  std::vector<SubMCEvent*>::size_type    siz = m_subMCEvents.size();
  if( 0 != siz ) {
    s << "\nSize of the Monte Carlo sub event vector :"
      << GlastEventField( GlastEvent::field4 )
      << siz;
    long count = 0;
    std::vector<SubMCEvent*>::const_iterator iter;
    for( iter = m_subMCEvents.begin(), count = 0;
              iter!=m_subMCEvents.end();
                    iter++, count++ ) {
      s << "\nIndex "
        << GlastEventField( GlastEvent::field4 )
        << count
        << " of object of type " << **iter;
    }
  }
  return s;
}


#endif    // GLASTEVENT_MCEVENT_H
