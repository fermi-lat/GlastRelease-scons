// $Header$
#ifndef LHCBEVENT_RANDOMNUMBERSEED_H
#define LHCBEVENT_RANDOMNUMBERSEED_H 1


// Include files
#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/Definitions.h"


//------------------------------------------------------------------------------
//
// ClassName:   RandomNumberSeed
//  
// Description: Random number seed for Monte Carlo
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class RandomNumberSeed                                                         {

public:

  /// Constructors
  RandomNumberSeed()
    : m_seed1(0),
      m_seed2(0)                                                             { }
  RandomNumberSeed (unsigned long seed1, unsigned long seed2)
    : m_seed1(seed1),
      m_seed2(seed2)                                                         { }
  /// Destructor
  ~RandomNumberSeed()                                                        { }

  /// Equality operator
  bool operator == (const RandomNumberSeed& t) const                           {
    return (m_seed1 == t.m_seed1) && (m_seed2 == t.m_seed2);
  }
  /// Retrieve m_seed1
  unsigned long seed1() const                                                  {
    return m_seed1;
  }
  /// Update m_seed1
  void setSeed1( unsigned long value )                                         {
    m_seed1 = value;
  }
  /// Retrieve m_seed2
  unsigned long seed2() const                                                  {
    return m_seed2;
  }
  /// Update m_seed2
  void setSeed2( unsigned long value )                                         {
    m_seed2 = value;
  }

  /// Serialize the object for writing
  friend StreamBuffer& operator<< ( StreamBuffer& s,
                                    const RandomNumberSeed& obj )              {
    return s << obj.m_seed1 << obj.m_seed2;
  }
  /// Serialize the object for reading
  friend StreamBuffer& operator>> ( StreamBuffer& s,
                                    RandomNumberSeed& obj )                    {
    return s >> obj.m_seed1 >> obj.m_seed2;
  }

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s,
                                    const RandomNumberSeed& obj )              {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  std::ostream& fillStream( std::ostream& s ) const                            {
    return s   << "class RandomNumberSeed ( seed1, seed2 ) : ( "
      << GlastEventField( GlastEvent::field12 )
      << m_seed1 << ", "
      << GlastEventField( GlastEvent::field12 )
      << m_seed2 << " )";
  }

private:

  /// Seed1
  unsigned long m_seed1;
  /// Seed2
  unsigned long m_seed2;

};


#endif    // LHCBEVENT_RANDOMNUMBERSEED_H
