#ifndef ROOT_TIMETONE_H
#define ROOT_TIMETONE_H 1

#include <iostream>

#include "TObject.h"
#include "GemTime.h"

/** @class TimeTone
* @brief encapsulate the TimeTone Markers that are put in the science data stream.
*
*       
*
* $Header$
*/

class TimeTone : public TObject {
  
public:

  /// Default c'tor.  Assigns sentinel values to all fields
  TimeTone()
    :m_incomplete(LSF_INVALID_UINT), m_timeSecs(LSF_INVALID_UINT),
     m_flywheeling(LSF_INVALID_UINT), m_flags(LSF_INVALID_UCHAR), m_timeHack() {
  }

  /// Standard c'tor.  Takes input values for all fields
  TimeTone( UInt_t incomplete, UInt_t timeSecs,
	    UInt_t flywheeling, UChar_t flags,
	    const GemTime& timeHack )
    :m_incomplete(incomplete), m_timeSecs(timeSecs), 
     m_flywheeling(flywheeling), m_flags(flags), m_timeHack(timeHack) {
  }

  /// Copy c'tor.  Nothing fancy, just copy all values.
  TimeTone(const TimeTone& other)
    :TObject(other), 
     m_incomplete(other.incomplete()), m_timeSecs(timeSecs()), 
     m_flywheeling(other.flywheeling()), m_flags(flags()), m_timeHack(other.timeHack()){
  }

  /// D'tor.  Nothing special.
  virtual ~TimeTone() {}
  
  /// Assignement operator.  Nothing fancy, just copy all values.
  inline TimeTone& operator=( const TimeTone& other ) {
    initialize(other.incomplete(),other.timeSecs(),
	       other.flywheeling(),other.flags(),
	       other.timeHack());
    return *this;
  }
  
  /// If this is non-zero part of the time tone is missing, check status bits 
  inline UInt_t incomplete() const { return m_incomplete; }
  
  /// Number of seconds since Epoch start at most recent time hack
  inline UInt_t timeSecs() const { return m_timeSecs; }
  
  /// Number of time tones since last complete time tone
  inline UInt_t flywheeling() const { return m_flywheeling; }
  
  /// All the flags at once
  inline UChar_t flags() const { return m_flags; }
  
  /// Make sure that the flags are valid
  inline Bool_t flagsValid() const { return m_flags != LSF_INVALID_UCHAR; }

  /// NO GPS lock, message w.r.t. LAT clock
  inline Bool_t missingGps() const { return (m_flags & enums::Lsf::TimeTone::MISSING_GPS_MASK) != 0; }
  
  /// NO 1-PPS signal at CPU level
  inline Bool_t missingCpuPps() const { return (m_flags & enums::Lsf::TimeTone::MISSING_CPU_MASK) != 0; }
  
  /// NO 1-PPS signal at LAT level
  inline Bool_t missingLatPps() const {  return (m_flags & enums::Lsf::TimeTone::MISSING_LAT_MASK) != 0; }
  
  /// NO 1-PPS signal at Spacecraft 
  inline Bool_t missingTimeTone() const  { return (m_flags & enums::Lsf::TimeTone::MISSING_TIMETONE_MASK) != 0; }

  /// The value of the GemTime registers at the time tone
  inline const GemTime& timeHack() const { return m_timeHack; }
  
  /// set everything at once
  inline void initialize( UInt_t incomplete, UInt_t timeSecs,
			  UInt_t flywheeling, UChar_t flags,
			  const GemTime& timeHack ) {
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
    if ( m_flags == LSF_INVALID_UCHAR ) m_flags = 0;
    if ( value ) { m_flags |= enums::Lsf::TimeTone::MISSING_GPS_MASK; }
    else { m_flags &= ~(enums::Lsf::TimeTone::MISSING_GPS_MASK) ; }
  }
  inline void setMissingCpu( Bool_t value ) {
    if ( m_flags == LSF_INVALID_UCHAR ) m_flags = 0;
    if ( value ) { m_flags |= enums::Lsf::TimeTone::MISSING_CPU_MASK; }
    else { m_flags &= ~(enums::Lsf::TimeTone::MISSING_CPU_MASK) ; }
  }
  inline void setMissingLAT( Bool_t value ) {
    if ( m_flags == LSF_INVALID_UCHAR ) m_flags = 0;
    if ( value ) { m_flags |= enums::Lsf::TimeTone::MISSING_LAT_MASK; }
    else { m_flags &= ~(enums::Lsf::TimeTone::MISSING_LAT_MASK); }
  }
  inline void setMissingTimeTone( Bool_t value ) {
    if ( m_flags == LSF_INVALID_UCHAR ) m_flags = 0;
    if ( value ) { m_flags |= enums::Lsf::TimeTone::MISSING_TIMETONE_MASK; }
    else { m_flags &= ~(enums::Lsf::TimeTone::MISSING_TIMETONE_MASK); }
  }

  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_incomplete = LSF_INVALID_UINT;
    m_timeSecs = LSF_INVALID_UINT;
    m_flywheeling = LSF_INVALID_UINT;
    m_flags = LSF_INVALID_UCHAR;
    m_timeHack.Clear("");
  }
  /// ROOT print function
  void Print(Option_t* /* option="" */) const;

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const TimeTone &, const std::string & name = "" ) const ; // for tests

private:
  
  /// If this is non-zero part of the time tone is missing, check status bits 
  UInt_t m_incomplete;
  /// Number of seconds since Epoch start at most recent time hack
  UInt_t m_timeSecs;
  /// Number of time tones since last complete time tone
  UInt_t m_flywheeling;
  /// missing signals [ GPS | 1-pps (CPU) | 1-pps (LAT) | 1-pps (SC) ]
  UChar_t m_flags;

  /// The value of the GemTime registers at the time tone
  GemTime m_timeHack;          

  ClassDef(TimeTone,1) //   Information about the 1-PPS TimeTone and GPS signal status
  
};

#endif    // ROOT_TIMETONE_H
