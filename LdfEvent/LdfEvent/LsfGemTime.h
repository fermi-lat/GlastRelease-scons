#ifndef LSF_GEMTIME_H
#define LSF_GEMTIME_H 1

#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"

/** @class Timehack
* @brief FIXME
*
* $Header$
*/

namespace LsfEvent {
 
  class GemTime {
    
  public:
    
    GemTime()
      : m_hacks(0), m_ticks(0){}
    
    GemTime( unsigned int hacks, unsigned int ticks )
      : m_hacks(hacks), m_ticks(ticks){}

    GemTime( const GemTime& other )
      : m_hacks(other.hacks()), m_ticks(other.ticks()){}
    
    ~GemTime() {
    }

    /// Assignement operator
    inline GemTime& operator=( const GemTime& other ) {
      set(other.hacks(),other.ticks());
      return *this;
    }
    
    /// number of time hacks (1-second ticks) at sample time
    inline unsigned int hacks() const { return m_hacks; }
    
    /// number of time ticks (50ns second ticks) at sample time
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
     
    /// Number of 1-PPS time hacks
    unsigned int m_hacks;
       
    /// Number of ticks of 20MhZ clock
    unsigned int m_ticks;
  }; 

}

#endif    // EVENT_TIMEHACK_H
