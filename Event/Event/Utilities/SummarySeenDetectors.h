// $Header$
#ifndef LHCBEVENT_SUMMARYSEENDETECTORS_H
#define LHCBEVENT_SUMMARYSEENDETECTORS_H 1


// Include files
#include <iostream>
#include "Gaudi/Kernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/Definitions.h"


//------------------------------------------------------------------------------
//
// ClassName:   SummarySeenDetectors
//  
// Description: Class containg summary of all detector elements
//              used to contruct the CandidateParticle
//              For the moment it is a placeholder containig the
//              "AXTK detector mask"
//
// Author:      G.Corti 
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class SummarySeenDetectors                                                     {

public:

  /// Constructors
  SummarySeenDetectors( long detectCode )
    : m_detectCode(detectCode)                                               { }
  SummarySeenDetectors()
    : m_detectCode(0)                                                        { }
  /// Destructor
  ~SummarySeenDetectors()                                                    { }

  /// Retrieve detector summary 
  long detectCode() const                                                      {
    return m_detectCode;
  }
  /// Update detector summary
  void setDetectCode(long value )                                              {
    m_detectCode = value;
  }

  /// Serialize the object for writing
  friend StreamBuffer& operator<< ( StreamBuffer& s,
                                    const SummarySeenDetectors& obj )          {
    return s << obj.m_detectCode;
  }
  /// Serialize the object for reading
  friend StreamBuffer& operator>> ( StreamBuffer& s,
                                    SummarySeenDetectors& obj )                {
    return s >> obj.m_detectCode;
  }

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s,
                                    const SummarySeenDetectors& obj )          {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  std::ostream& fillStream( std::ostream& s ) const                            {
    return s << "class SummarySeenDetectors : "
      << GlastEventField( GlastEvent::field5 )
      << m_detectCode;
  }

private:

  /// Coded Summary Detector word
  long m_detectCode;

};


#endif    // LHCBEVENT_SUMMARYSEENDETECTORS_H
