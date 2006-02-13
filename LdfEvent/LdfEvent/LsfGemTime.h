#ifndef LSF_GEMTIME_H
#define LSF_GEMTIME_H 1

#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"

#include "enums/Lsf.h"

/** @class Timehack
* @brief FIXME
*
* $Header$
*/

namespace LsfEvent {

  class GemTime {

  public:
    
    /// Default c'tor.  Assigns sentinel values to all fields  
    GemTime()
      : m_hacks(LSF_INVALID_UINT), m_ticks(LSF_INVALID_UINT){}
    
    /// Standard c'tor.  Takes input values for all fields
    GemTime( unsigned int hacks, unsigned int ticks )
      : m_hacks(hacks), m_ticks(ticks){}
    
    /// Copy c'tor.  Nothing fancy, just copy all values.
    GemTime( const GemTime& other )
      : m_hacks(other.hacks()), m_ticks(other.ticks()){}
    
    /// D'tor.  Nothing special.
    virtual ~GemTime() {
    }
    
    /// Assignement operator.  Nothing fancy, just copy all values.
    inline GemTime& operator=( const GemTime& other ) {
      set(other.hacks(),other.ticks());
      return *this;
    }
    
    /// number of 1-PPS time hacks.  This counter wraps around at 127)
    inline unsigned int hacks() const { return m_hacks; }
    
  /// value of the LAT system clock at the last 1-PPS time hack
    inline unsigned int ticks() const { return m_ticks; }
    
    /// set everything at once
    inline void set(unsigned int hacks, unsigned int ticks) {
      m_hacks = hacks;
      m_ticks = ticks;
    }
    
    /// Update time ticks
    inline void setTicks( unsigned int value ) { m_ticks = value; }
    
    /// Update time hacks
    inline void setHacks( unsigned int value ) { m_hacks = value; }
        
    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const GemTime& obj )    {
      return s << obj.m_hacks << ' ' << obj.m_ticks;
    }
    /// Serialize the object for reading
    friend StreamBuffer& operator>> ( StreamBuffer& s, GemTime& obj )          {
      return s >> obj.m_hacks >> obj.m_ticks;
    }

    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const GemTime& obj )    {
      return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const                            {
      return s << "class GemTime : "
	       << EventField( EventFormat::field12 )
	       << m_hacks << ':' << m_ticks;
    }
    
  private:
     
    /// Number of 1-PPS time hacks, this is a simple 7-bit counter which rolls-over
    unsigned int m_hacks;
       
    /// Number of ticks of 20MhZ clock at last 1-PPS time hack
    unsigned int m_ticks;
  }; 

}

#endif    // EVENT_TIMEHACK_H
