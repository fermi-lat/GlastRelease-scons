#ifndef LSFDATA_TIME_H
#define LSFDATA_TIME_H 1

#include <iostream>
#include <string>
#include <cstdio>

#include "lsfData/LsfTimeTone.h"
#include "lsfData/LsfGemTime.h"

/** @class Time
*
* $Header$
*/

namespace lsfData {

  class Time {
    
  public:
    
    Time( const TimeTone& current, const TimeTone& previous,
	  const GemTime& timeHack, unsigned int timeTicks) 
      :m_current(current),m_previous(previous),
       m_timeHack(timeHack),m_timeTicks(timeTicks) {
    }
    
    Time()
      :m_current(),m_previous(),
       m_timeHack(),m_timeTicks(0){
    }
    
    Time( const Time& other )
      :m_current(other.current()),m_previous(other.previous()),
       m_timeHack(other.timeHack()),m_timeTicks(other.timeTicks()){
    }
    
    virtual ~Time(){
    }
    
    /// Assignement operator
    inline Time& operator=( const Time& other ) {
      set(other.current(),other.previous(),
	  other.timeHack(),other.timeTicks());
      return *this;
    }

    void clear() {
       m_current.clear();
       m_previous.clear();
       m_timeHack.clear();
       m_timeTicks = 0;
    }
   
    void print() const {
      printf(" trigger:  tics = 0x%08X (%u)\n", m_timeTicks, m_timeTicks );
      m_timeHack.print(" ");
        printf(" current:  secs = %u\n", m_current.timeSecs());
        m_current.timeHack().print(" ");
        printf(" previous:  secs = %u\n", m_previous.timeSecs());
        m_previous.timeHack().print(" ");
    }

    friend std::ostream& operator<< ( std::ostream& s, const Time& obj )
    {
      return obj.fillStream(s);
    }

    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const {

      s << " current:  secs = " << m_current.timeSecs() << "\n"
        << "   tics  = 0x" << std::hex << std::uppercase << std::setfill('0') 
        << std::setw(8)
        << m_current.timeHack().ticks() << " (" 
        << std::dec << m_current.timeHack().ticks() << ")" << std::endl
        << "   hacks = 0x" << std::hex << m_current.timeHack().hacks() << " ("
        << std::dec << m_current.timeHack().hacks() << ")" << std::endl
        << " previous:  secs = " << m_previous.timeSecs() << "\n"
        << "   tics  = 0x" << std::hex << m_previous.timeHack().ticks() << " (" 
        << std::dec << m_previous.timeHack().ticks() << ")" << std::endl
        << "   hacks = 0x" << std::hex << m_previous.timeHack().hacks() << " ("
        << std::dec << m_previous.timeHack().hacks() << ")" << std::endl;

      return s;
    }

    void printHack(const std::string str="") const {
        printf("%stimeTics = 0x%08X (%u)\n", str.c_str(), m_timeTicks, m_timeTicks);
        m_timeHack.print(str);
    }

    /// The TimeTone right before the "active" one at event capture time
    inline const TimeTone& current() const { return m_current; } 

    /// The TimeTone that was "active" at event capture time
    inline const TimeTone& previous() const { return m_previous; }

    /// The GemTime at event capture time
    inline const GemTime& timeHack() const { return m_timeHack; } 

    /// The number of 50ns ticks since last the last time hack
    inline unsigned int timeTicks() const { return m_timeTicks; }

    /// set everything at once
    inline void set(const TimeTone& current, const TimeTone& previous,
		    const GemTime& timeHack, unsigned int timeTicks) {
      m_current = current;
      m_previous = previous;      
      m_timeHack = timeHack; 
      m_timeTicks = timeTicks;
    }

    // set the individual data members
    inline void setCurrent( const TimeTone& val) { m_current = val; }; 
    inline void setPrevious( const TimeTone& val) { m_previous = val; };
    inline void setGemTime( const GemTime& val) { m_timeHack = val; }; 
    inline void setTimeTicks( unsigned int val) { m_timeTicks = val; };
    
    
  private:
    
    /// 
    TimeTone m_current;
    TimeTone m_previous;
    GemTime m_timeHack;
    unsigned int m_timeTicks;
     
  };

}


#endif   
