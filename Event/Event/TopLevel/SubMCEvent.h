// $Header$
#ifndef LHCBEVENT_SUBMCEVENT_H
#define LHCBEVENT_SUBMCEVENT_H 1


// Include files
#include <iostream>
#include <vector>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GlastEvent/Utilities/RandomNumberSeed.h"
#include "GlastEvent/Utilities/ProcessingVersion.h"
#include "CLHEP/Geometry/Point3D.h"
#include "GlastEvent/Utilities/CLHEPStreams.h"
#include "GlastEvent/TopLevel/Definitions.h"


//------------------------------------------------------------------------------
//
// ClassName:   SubMCEvent
//  
// Description: Essential information of the sub- Monte Carlo event
//              (the Monte Carlo event can be assembled from several
//              of these sub-events)
//
// The class SubMCEvent uses the Class Library for HEP
// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/index.html">CLHEP</A>
//              
//              It contains:
//              - detector characteristics
//              - primary vertex
//              - Monte Carlo event weight
//              - random number seed
//              - processing version
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class SubMCEvent                                                               {

public:
  /// Constructors
  SubMCEvent()                                                               { }
  /// Destructor
  ~SubMCEvent()                                                              { }

  /// Clone operator
  SubMCEvent& operator= ( const SubMCEvent clone )                             {
    m_detectorCharacteristics = clone.detectorCharacteristics();
    m_processingVersion = clone.processingVersion();
    m_randomNumberSeed = clone.randomNumberSeed();
    m_primaryVertex = clone.primaryVertex();
    m_weight = clone.weight();
    return *this;
  }

  /// Equality operator
  bool operator== ( const SubMCEvent& t ) const                                {
    return 
      detectorCharacteristics() == t.detectorCharacteristics()
      && processingVersion() == t.processingVersion()
      && randomNumberSeed() == t.randomNumberSeed()
      && primaryVertex() == t.primaryVertex()
      && weight() == t.weight();
  }

  /// Retrieve detectorCharacteristics
  long detectorCharacteristics () const                                        {
    return m_detectorCharacteristics;
  }
  /// Update detectorCharacteristics
  void setDetectorCharacteristics (long value)                                 {
    m_detectorCharacteristics = value;
  }

  /// Retrieve primary vertex location
  const HepPoint3D& primaryVertex () const                                     {
    return m_primaryVertex;
  }
  HepPoint3D& primaryVertex ()                                                 {
    return m_primaryVertex;
  }
  /// Update primary vertex location
  void setPrimaryVertex (const HepPoint3D& value)                              {
    m_primaryVertex = value;
  }

  /// Retrieve Monte Carlo event weight
  double weight () const                                                       {
    return m_weight;
  }
  /// Update Monte Carlo event weight
  void setWeight (double value)                                                {
    m_weight = value;
  }

  /// Retrieve Monte Carlo random number generator seed
  const RandomNumberSeed& randomNumberSeed () const                            {
    return m_randomNumberSeed; 
  }
  /// Update Monte Carlo random number generator seed
  void setRandomNumberSeed (const RandomNumberSeed& value)                     {
    m_randomNumberSeed = value;
  }

  /// Retrieve processing version
  const ProcessingVersion& processingVersion () const                          {
    return m_processingVersion; 
  }
  /// Update processing version
  void setProcessingVersion (const ProcessingVersion& value)                   {
    m_processingVersion = value; 
  }

  /// Serialize the object for writing
  friend StreamBuffer& operator<< ( StreamBuffer& s, const SubMCEvent& obj )   {
    return s
      << obj.m_detectorCharacteristics
      << obj.m_primaryVertex
      << obj.m_weight
      << obj.m_randomNumberSeed
      << obj.m_processingVersion;
  }
  /// Serialize the object for reading
  friend StreamBuffer& operator>> ( StreamBuffer& s, SubMCEvent& obj )         {
    return s
      >> obj.m_detectorCharacteristics
      >> obj.m_primaryVertex
      >> obj.m_weight
      >> obj.m_randomNumberSeed
      >> obj.m_processingVersion;
  }


  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const SubMCEvent& obj )   {
    return obj.fillStream(s);
  }  
  /// Fill the output stream (ASCII)
  std::ostream& fillStream( std::ostream& s ) const;

private: 
  /// Detector characteristics
  long                m_detectorCharacteristics;
  /// Primary vertex location
  /// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/Geometry/HepPoint3D.html">class HepPoint3D</A>
  HepPoint3D          m_primaryVertex;
  /// Monte Carlo event weight
  double              m_weight;
  /// Monte Carlo random number generator seed
  RandomNumberSeed    m_randomNumberSeed;
  /// Processing version
  ProcessingVersion   m_processingVersion;

};


//
// Inline code must be outside the class definition
//


/// Fill the output stream (ASCII)
inline std::ostream& SubMCEvent::fillStream( std::ostream& s ) const           {
  return s
    << "class SubMCEvent :\n"
    << "\n    Detector characteristics = "
    << GlastEventField( GlastEvent::field12 )
    << m_detectorCharacteristics
    << "\n    Primary vertex location  = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_primaryVertex.x() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_primaryVertex.y() << ", "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_primaryVertex.z() << " )"
    << "\n    Monte Carlo event weight = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_weight
    << "\n    Random number seed       = " << m_randomNumberSeed
    << "\n    Processing version       = " << m_processingVersion;
}


#endif    // LHCBEVENT_SUBMCEVENT_H
