#ifndef ROOT_TIMETONE_H
#define ROOT_TIMETONE_H 1

#include <iostream>

#include "TObject.h"
#include "LsfTimeHack.h"

/** @class TimeTone
* @brief encapsulate the time tone markers that come down with the 
*
* $Header$
*/

class TimeTone : public TObject {
  
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
    
  TimeTone( UInt_t incomplete, UInt_t timeSecs,
	    UInt_t flywheeling, UChar_t flags,
	    const TimeHack& timeHack )
    :m_incomplete(incomplete), m_timeSecs(timeSecs), 
     m_flywheeling(flywheeling), m_flags(flags), m_timeHack(timeHack) {
  }
  
  TimeTone(const TimeTone& other)
    :TObject(other), 
     m_incomplete(other.incomplete()), m_timeSecs(timeSecs()), 
     m_flywheeling(other.flywheeling()), m_flags(flags()), m_timeHack(other.timeHack()){
  }
  
  ~TimeTone() {}
  
  /// Assignement operator
  inline TimeTone& operator=( const TimeTone& other ) {
    initialize(other.incomplete(),other.timeSecs(),
	       other.flywheeling(),other.flags(),
	       other.timeHack());
    return *this;
  }
  
  /// If this is non-zero part of the time tone is missing, check status bits 
  inline UInt_t incomplete() const { return m_incomplete; }
  
  /// Number of seconds since Epoch start at time hack
  inline UInt_t timeSecs() const { return m_timeSecs; }
  
  /// Number of time tones since last complete time tone
  inline UInt_t flywheeling() const { return m_flywheeling; }
  
  /// All the flags at once
  inline UChar_t flags() const { return m_flags; }
  
  /// NO GPS lock, message w.r.t. LAT clock
  inline Bool_t missingGps() const { return (m_flags & MISSING_GPS_MASK) != 0; }
  
  /// NO 1-PPS signal at CPU level
  inline Bool_t missingCpuPps() const { return (m_flags & MISSING_CPU_MASK) != 0; }
  
  /// NO 1-PPS signal at LAT level
  inline Bool_t missingLatPps() const {  return (m_flags & MISSING_LAT_MASK) != 0; }
  
  /// NO 1-PPS signal at Spacecraft 
  inline Bool_t missingTimeTone() const  { return (m_flags & MISSING_TIMETONE_MASK) != 0; }

  /// The time hack at the time tone
  inline const TimeHack& timeHack() const { return m_timeHack; }
  
  /// set everything at once
  inline void initialize( UInt_t incomplete, UInt_t timeSecs,
			  UInt_t flywheeling, UChar_t flags,
			  const TimeHack& timeHack ) {
    m_incomplete = incomplete;
    m_timeSecs = timeSecs;
    m_flywheeling = flywheeling;
    m_flags = flags;
    m_timeHack = timeHack;
  }
  
  // set the individual data members
  inline void setIncomplete( UInt_t value ) { m_incomplete = value; }
  inline void setTimeSecs( UInt_t value ) { m_timeSecs = value; }
  inline void setFlywheeling( UInt_t value ) { m_flywheeling = value; }  
  inline void setFlags( UChar_t value ) { m_flags = value; }
  inline void setMissingGps( Bool_t value ) { 
    if ( value ) { m_flags |= MISSING_GPS_MASK; }
    else { m_flags &= ~(MISSING_GPS_MASK) ; }
  }
  inline void setMissingCpu( Bool_t value ) {
    if ( value ) { m_flags |= MISSING_CPU_MASK; }
    else { m_flags &= ~(MISSING_CPU_MASK) ; }
  }
  inline void setMissingLAT( Bool_t value ) {
    if ( value ) { m_flags |= MISSING_LAT_MASK; }
    else { m_flags &= ~(MISSING_LAT_MASK); }
  }
  inline void setMissingTimeTone( Bool_t value ) {
    if ( value ) { m_flags |= MISSING_TIMETONE_MASK; }
    else { m_flags &= ~(MISSING_TIMETONE_MASK); }
  }

  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_incomplete = 0;
    m_timeSecs = 0;
    m_flywheeling = 0;
    m_flags = 0;
    m_timeHack.Clear("");
  }
  /// ROOT print function
  void Print(Option_t* /* option="" */) const;

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const TimeTone &, const std::string & name = "" ) const ; // for tests

private:
  
  /// Time
  UInt_t m_incomplete;  // 0 -> ok
  UInt_t m_timeSecs;    // # of secs since epoch @ time hack
  UInt_t m_flywheeling; // # of timetones since last complete once
  unsigned char m_flags;       // missing signals [ GPS | 1-pps (CPU) | 1-pps (LAT) | 1-pps (SC) ]
  
  TimeHack m_timeHack;          

  ClassDef(TimeTone,1) ;
  
};

#endif    // ROOT_TIMETONE_H
