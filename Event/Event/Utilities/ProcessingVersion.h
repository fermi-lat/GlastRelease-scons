// $Header$
#ifndef LHCBEVENT_PROCESSINGVERSION_H
#define LHCBEVENT_PROCESSINGVERSION_H 1


// Include files
#include <iostream>
#include "Gaudi/Kernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/Definitions.h"


//------------------------------------------------------------------------------
//
// ClassName:   ProcessingVersion
//  
// Description: Placeholder for processing version (used on the (sub)event level
//              and on the level of LHCb containers of physics objects)
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class ProcessingVersion                                                        {

public:

  /// Constructors
  ProcessingVersion( long version )
    : m_version(version)                                                     { }
  ProcessingVersion()
    : m_version(0)                                                           { }
  /// Destructor
  ~ProcessingVersion()                                                       { }

  /// Retrieve version number
  long versionNumber() const                                                   {
    return m_version;
  }
  /// Update version number
  void setVersionNumber( long value )                                          {
    m_version = value;
  }
  /// Ordering operator
  friend bool operator < ( const ProcessingVersion& left,
                           const ProcessingVersion& right)                     {
    return left.m_version < right.m_version;
  }
  /// Equality operator
  friend bool operator == ( const ProcessingVersion& left,
                            const ProcessingVersion& right)                    {
    return left.m_version == right.m_version;
  }

  /// Serialize the object for writing
  friend StreamBuffer& operator<< ( StreamBuffer& s,
                                    const ProcessingVersion& obj )             {
    return s << obj.m_version;
  }
  /// Serialize the object for reading
  friend StreamBuffer& operator>> ( StreamBuffer& s,
                                    ProcessingVersion& obj )                   {
    return s >> obj.m_version;
  }

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s,
                                    const ProcessingVersion& obj )             {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  std::ostream& fillStream( std::ostream& s ) const                            {
    return s << "class ProcessingVersion : "
      << GlastEventField( GlastEvent::field5 )
      << m_version;
  }

private:

  // Version number
  long m_version;

};


#endif    // LHCBEVENT_PROCESSINGVERSION_H
