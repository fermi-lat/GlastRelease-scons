// $Header$
#ifndef LHCBEVENT_DETECTORDATAOBJECT_H
#define LHCBEVENT_DETECTORDATAOBJECT_H 1


// Include files
#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/Definitions.h"


//------------------------------------------------------------------------------
//
// ClassName:   DetectorDataObject
//  
// Description: Placeholder for smart pointer to the detector description
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class DetectorDataObject                                                       {

public:

  /// Constructors
  DetectorDataObject (long detector)
    : m_detector(detector)                                                   { }
  DetectorDataObject()
    : m_detector(0)                                                          { }
  /// Destructor
  ~DetectorDataObject()                                                      { }

  /// Retrieve detector data object
  long detector() const                                                        {
    return m_detector;
  }
  /// Update detector data object
  void setDetector( long value )                                               {
    m_detector = value;
  }

  /// Serialize the object for writing
  friend StreamBuffer& operator<< ( StreamBuffer& s,
                                    const DetectorDataObject& obj )            {
    return s << obj.m_detector;
  }
  /// Serialize the object for reading
  friend StreamBuffer& operator>> ( StreamBuffer& s,
                                    DetectorDataObject& obj )                  {
    return s >> obj.m_detector;
  }

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s,
                                    const DetectorDataObject& obj )            {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  std::ostream& fillStream( std::ostream& s ) const                            {
    return s << "class DetectorDataObject : "
      << GlastEventField( GlastEvent::field5 )
      << m_detector;
  }

private:

  // Detector data object identifier
  long m_detector;

};


#endif  // LHCBEVENT_DETECTORDATAOBJECT_H
