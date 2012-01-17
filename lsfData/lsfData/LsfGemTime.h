#ifndef LSFDATA_GEMTIME_H
#define LSFDATA_GEMTIME_H 1

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>

/** @class Timehack
* @brief FIXME
*
* $Header$
*/

namespace lsfData {
 
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
    
    void clear() {
        m_hacks = 0;
        m_ticks = 0;
    }

    void print(std::string str="") const {
        printf("%stics  = 0x%08X (%u)\n", str.c_str(), m_ticks, m_ticks);
        printf("%shacks = 0x%08X (%u)\n", str.c_str(), m_hacks, m_hacks);
    }

    friend std::ostream& operator<< ( std::ostream& s, const GemTime& obj )
    {
      return obj.fillStream(s);
    }

    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const {
      
      s << "tics  = 0x" << std::hex << std::setfill('0') << std::setw(8)
        << m_ticks 
        << " (" << std::dec << m_ticks << ")\n"
        << "hacks = 0x" << std::hex << m_hacks << " (" << std::dec
        << m_hacks << ")\n";
      return s;
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
    
    
  private:
     
    /// Number of 1-PPS time hacks
    unsigned int m_hacks;
       
    /// Number of ticks of 20MhZ clock
    unsigned int m_ticks;
  }; 

}

#endif    // EVENT_TIMEHACK_H
