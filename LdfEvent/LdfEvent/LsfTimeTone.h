#ifndef LSF_TIMETONE_H
#define LSF_TIMETONE_H 1

#include <iostream>

#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"
#include "LdfEvent/LsfGemTime.h"

#include "lsfDataStore/LsfTimeTone.h"

/** @class TimeTone
* @brief FIXME
*
* $Header$
*/

namespace LsfEvent {
  
  class TimeTone {
    
  public:
    
    /// Default c'tor.  Assigns sentinel values to all fields
    TimeTone()
      :m_incomplete(LSF_INVALID_UINT), m_timeSecs(LSF_INVALID_UINT),
       m_flywheeling(LSF_INVALID_UINT), m_flags(LSF_INVALID_UCHAR), m_timeHack() {
    }
    
    /// Standard c'tor.  Takes input values for all fields
    TimeTone( unsigned int incomplete, unsigned int timeSecs,
	      unsigned int flywheeling, unsigned char flags,
	      const GemTime& timeHack )
      :m_incomplete(incomplete), m_timeSecs(timeSecs), 
       m_flywheeling(flywheeling), m_flags(flags), m_timeHack(timeHack) {
    }
    
    /// Copy c'tor.  Nothing fancy, just copy all values.
    TimeTone(const TimeTone& other)
      :m_incomplete(other.incomplete()), m_timeSecs(timeSecs()), 
       m_flywheeling(other.flywheeling()), m_flags(flags()), m_timeHack(other.timeHack()){
    }
    
    TimeTone(const lsfDataStore::TimeTone& other)
      :m_incomplete(other.incomplete()), m_timeSecs(other.timeSecs()), 
       m_flywheeling(other.flywheeling()), m_flags(other.flags()), m_timeHack(other.timeHack()){
    }
    
    /// D'tor.  Nothing special.
    virtual ~TimeTone() {}
    
    /// Assignement operator.  Nothing fancy, just copy all values.
    inline TimeTone& operator=( const TimeTone& other ) {
      set(other.incomplete(),other.timeSecs(),
	  other.flywheeling(),other.flags(),
	  other.timeHack());
      return *this;
    }
    
    /// If this is non-zero part of the time tone is missing, check status bits 
    inline unsigned int incomplete() const { return m_incomplete; }
    
    /// Number of seconds since Epoch start at most recent time hack
    inline unsigned int timeSecs() const { return m_timeSecs; }
    
    /// Number of time tones since last complete time tone
    inline unsigned int flywheeling() const { return m_flywheeling; }
    
    /// All the flags at once
    inline unsigned char flags() const { return m_flags; }
    
    /// Make sure that the flags are valid
    inline bool flagsValid() const { return m_flags != LSF_INVALID_UCHAR; }

    /// NO GPS lock, message w.r.t. LAT clock
    inline bool missingGps() const { return (m_flags & enums::Lsf::TimeTone::MISSING_GPS_MASK) != 0; }
  
    /// NO 1-PPS signal at CPU level
    inline bool missingCpuPps() const { return (m_flags & enums::Lsf::TimeTone::MISSING_CPU_MASK) != 0; }
    
  /// NO 1-PPS signal at LAT level
    inline bool missingLatPps() const {  return (m_flags & enums::Lsf::TimeTone::MISSING_LAT_MASK) != 0; }
    
    /// NO 1-PPS signal at Spacecraft 
    inline bool missingTimeTone() const  { return (m_flags & enums::Lsf::TimeTone::MISSING_TIMETONE_MASK) != 0; }
    
    /// The value of the GemTime registers at the time tone
    inline const GemTime& timeHack() const { return m_timeHack; }
    
    /// set everything at once
    inline void set( unsigned int incomplete, unsigned int timeSecs,
		     unsigned int flywheeling, unsigned char flags,
		     const GemTime& timeHack ) {
      m_incomplete = incomplete;
      m_timeSecs = timeSecs;
      m_flywheeling = flywheeling;
      m_flags = flags;
      m_timeHack = timeHack;
    }
    
    // set the individual data members
    inline void setIncomplete( unsigned int value ) { m_incomplete = value; }
    inline void setTimeSecs( unsigned int value ) { m_timeSecs = value; }
    inline void setFlywheeling( unsigned int value ) { m_flywheeling = value; }  
    inline void setFlags( unsigned char value ) { m_flags = value; }
    inline void setMissingGps( bool value ) { 
      if ( m_flags == LSF_INVALID_UCHAR ) m_flags = 0;
      if ( value ) { m_flags |= enums::Lsf::TimeTone::MISSING_GPS_MASK; }
      else { m_flags &= ~(enums::Lsf::TimeTone::MISSING_GPS_MASK) ; }
    }
    inline void setMissingCpu( bool value ) {
      if ( m_flags == LSF_INVALID_UCHAR ) m_flags = 0;
      if ( value ) { m_flags |= enums::Lsf::TimeTone::MISSING_CPU_MASK; }
      else { m_flags &= ~(enums::Lsf::TimeTone::MISSING_CPU_MASK) ; }
    }
    inline void setMissingLAT( bool value ) {
      if ( m_flags == LSF_INVALID_UCHAR ) m_flags = 0;
      if ( value ) { m_flags |= enums::Lsf::TimeTone::MISSING_LAT_MASK; }
      else { m_flags &= ~(enums::Lsf::TimeTone::MISSING_LAT_MASK); }
    }
    inline void setMissingTimeTone( bool value ) {
      if ( m_flags == LSF_INVALID_UCHAR ) m_flags = 0;
      if ( value ) { m_flags |= enums::Lsf::TimeTone::MISSING_TIMETONE_MASK; }
      else { m_flags &= ~(enums::Lsf::TimeTone::MISSING_TIMETONE_MASK); }
    }
    
    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const TimeTone& obj )    {
      s << obj.m_incomplete << ' ' << obj.m_timeSecs << ' ' << obj.m_flywheeling << ' ' << obj.m_flags;
      return s;
    }
    /// Serialize the object for reading
    friend StreamBuffer& operator>> ( StreamBuffer& s, TimeTone& obj )          {
      s >> obj.m_incomplete >> obj.m_timeSecs >> obj.m_flywheeling >> obj.m_flags;
      return s;
    }
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const TimeTone& obj )    {
      return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const                            {
      s << "class TimeTone : " 
	<< "Secs " << m_timeSecs << " Incomplete: " << m_incomplete 
        << "Flywheeling: "  << m_flywheeling << "Flags: " << m_flags;
      return s;
    }

  private:
    
    /// If this is non-zero part of the time tone is missing, check status bits 
    unsigned int m_incomplete;
    /// Number of seconds since Epoch start at most recent time hack
    unsigned int m_timeSecs;
    /// Number of time tones since last complete time tone
    unsigned int m_flywheeling;
    /// missing signals [ GPS | 1-pps (CPU) | 1-pps (LAT) | 1-pps (SC) ]
    unsigned char m_flags;
    
    /// The value of the GemTime registers at the time tone
    GemTime m_timeHack;          
    
  };

}

#endif    // EVENT_TIMETONE_H
