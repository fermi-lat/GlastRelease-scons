#ifndef LSF_TIME_H
#define LSF_TIME_H 1

#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"

#include "LdfEvent/LsfTimeTone.h"
#include "LdfEvent/LsfTimeHack.h"

/** @class Time
*
* $Header$
*/

namespace LsfEvent {

  class Time {
    
  public:
    
    Time( const TimeTone& current, const TimeTone& previous,
	  const TimeHack& timeHack, unsigned int timeTicks)
      :m_current(current),m_previous(previous),
       m_timeHack(timeHack),m_timeTicks(timeTicks){
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

    /// The TimeTone right before the "active" one at event capture time
    inline const TimeTone& current() const { return m_current; } 

    /// The TimeTone that was "active" at event capture time
    inline const TimeTone& previous() const { return m_previous; }

    /// The TimeHack at event capture time
    inline const TimeHack& timeHack() const { return m_timeHack; } 

    /// The number of 50ns ticks since last the last time hack
    inline unsigned int timeTicks() const { return m_timeTicks; }

    /// set everything at once
    inline void set(const TimeTone& current, const TimeTone& previous,
		    const TimeHack& timeHack, unsigned int timeTicks) {
      m_current = current;
      m_previous = previous;      
      m_timeHack = timeHack; 
      m_timeTicks = timeTicks;
    }

    // set the individual data members
    inline void setCurrent( const TimeTone& val) { m_current = val; }; 
    inline void setPrevious( const TimeTone& val) { m_previous = val; };
    inline void setTimeHack( const TimeHack& val) { m_timeHack = val; }; 
    inline void setTimeTicks( unsigned int val) { m_timeTicks = val; };
    
    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const Time& obj )    {
      s << obj.m_current; 
      s << obj.m_previous;
      s << obj.m_timeHack;
      s << obj.m_timeTicks;
      return s;
    }
    /// Serialize the object for reading
    friend StreamBuffer& operator>> ( StreamBuffer& s, Time& obj )          {
      s >> obj.m_current;
      s >> obj.m_previous;
      s >> obj.m_timeHack;
      s >> obj.m_timeTicks;
      return s;
    }
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const Time& obj )    {
      return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const                            {
      return s << "class Time : " << std::endl
	       << EventField( EventFormat::field12 )
	       << m_current << std::endl
	       << m_previous << std::endl
	       << m_timeHack << std::endl
	       << m_timeTicks;
    }
      
    
  private:
    
    /// 
    TimeTone m_current;
    TimeTone m_previous;
    TimeHack m_timeHack;
    unsigned int m_timeTicks;
     
  };

}


#endif    // EVENT_TIMETONE_H
