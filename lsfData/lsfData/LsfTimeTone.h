#ifndef LSFDATA_TIMETONE_H
#define LSFDATA_TIMETONE_H 1

#include <iostream>

#include "lsfData/LsfGemTime.h"

/** @class TimeTone
* @brief FIXME
*
* $Header$
*/

namespace lsfData {
  
  class TimeTone {
    
  public:
    
    enum { MISSING_GPS_MASK = 1, 
	   MISSING_CPU_MASK = 2,
	   MISSING_LAT_MASK = 4,
	   MISSING_TIMETONE_MASK = 8 } MASKS;
    
  public:
    
    TimeTone()
      :m_incomplete(0), m_timeSecs(0), 
       m_flywheeling(0), m_flags(0), m_timeHack() {
    }
    
    TimeTone( unsigned int incomplete, unsigned int timeSecs,
	      unsigned int flywheeling, unsigned char flags,
	      const GemTime& timeHack )
    :m_incomplete(incomplete), m_timeSecs(timeSecs), 
     m_flywheeling(flywheeling), m_flags(flags), m_timeHack(timeHack) {
    }
    
    TimeTone(const TimeTone& other)
      :m_incomplete(other.incomplete()), m_timeSecs(other.timeSecs()), 
       m_flywheeling(other.flywheeling()), m_flags(other.flags()), m_timeHack(other.timeHack()){
    }

    ~TimeTone() {}
    
    /// Assignement operator
    inline TimeTone& operator=( const TimeTone& other ) {
      set(other.incomplete(),other.timeSecs(),
	  other.flywheeling(),other.flags(),
	  other.timeHack());
      return *this;
    }
    
    void clear() {
        m_incomplete = 0;
        m_timeSecs = 0;
        m_flywheeling = 0;
        m_flags = 0;
        m_timeHack.clear();
    }

    /// If this is non-zero part of the time tone is missing, check status bits 
    inline unsigned int incomplete() const { return m_incomplete; }
    
    /// Number of seconds since Epoch start at time hack
    inline unsigned int timeSecs() const { return m_timeSecs; }
    
    /// Number of time tones since last complete time tone
    inline unsigned int flywheeling() const { return m_flywheeling; }
    
    /// All the flags at once
    inline unsigned char flags() const { return m_flags; }
    
    /// NO GPS lock, message w.r.t. LAT clock
    inline bool missingGps() const { return (m_flags & MISSING_GPS_MASK) != 0; }
    
    /// NO 1-PPS signal at CPU level
    inline bool missingCpuPps() const { return (m_flags & MISSING_CPU_MASK) != 0; }
    
    /// NO 1-PPS signal at LAT level
    inline bool missingLatPps() const {  return (m_flags & MISSING_LAT_MASK) != 0; }
    
    /// NO 1-PPS signal at Spacecraft 
    inline bool missingTimeTone() const  { return (m_flags & MISSING_TIMETONE_MASK) != 0; }

    /// The time hack
    inline const GemTime& timeHack() const { return m_timeHack; }
    
    /// set everything at once
    inline void set( unsigned int incomplete, unsigned int timeSecs,
		     unsigned int flywheeling, unsigned char flags,
		     const GemTime& timeHack ) {
      m_incomplete = incomplete;
      m_timeSecs = timeSecs;
      m_flywheeling = flywheeling;
      m_flags = flags;
      m_timeHack.set(timeHack.hacks(),timeHack.ticks());
    }
    
    // set the individual data members
    inline void setIncomplete( unsigned int value ) { m_incomplete = value; }
    inline void setTimeSecs( unsigned int value ) { m_timeSecs = value; }
    inline void setFlywheeling( unsigned int value ) { m_flywheeling = value; }  
    inline void setFlags( unsigned char value ) { m_flags = value; }
    inline void setMissingGps( bool value ) { 
      if ( value ) m_flags |= MISSING_GPS_MASK; 
      else m_flags &= (~MISSING_GPS_MASK);
    }
    inline void setMissingCpu( bool value ) {
      if ( value ) m_flags |= MISSING_CPU_MASK;
      else m_flags &= (~MISSING_CPU_MASK);
    }
    inline void setMissingLAT( bool value ) {
      if ( value ) m_flags |= MISSING_LAT_MASK;
      else m_flags &= (~MISSING_LAT_MASK);
    }
    inline void setMissingTimeTone( bool value ) {
      if ( value ) m_flags |= MISSING_TIMETONE_MASK;
      else m_flags &= (~MISSING_TIMETONE_MASK);
    }


  private:
    
    /// Time
    unsigned int m_incomplete;  // 0 -> ok
    unsigned int m_timeSecs;    // # of secs since epoch @ time hack
    unsigned int m_flywheeling; // # of timetones since last complete once
    unsigned char m_flags;       // missing signals [ GPS | 1-pps (CPU) | 1-pps (LAT) | 1-pps (SC) ]
    
    GemTime m_timeHack;          
    
  };

}

#endif    // EVENT_TIMETONE_H
