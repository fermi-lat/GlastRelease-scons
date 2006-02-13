#ifndef ROOT_TIME_H
#define ROOT_TIME_H 1

#include <iostream>

#include "LsfTimeTone.h"
#include "LsfGemTime.h"

/** @class Time
*
* $Header$
*/


class Time : public TObject {
  
public:

  /// Default c'tor.  Assigns sentinel values to all fields
  Time()
    :m_current(),m_previous(),
    m_timeHack(),m_timeTicks(LSF_INVALID_UINT){
  }
  
  /// Standard c'tor.  Takes input values for all fields
  Time( const TimeTone& current, const TimeTone& previous,
	const GemTime& timeHack, UInt_t timeTicks)
    :m_current(current),m_previous(previous),
     m_timeHack(timeHack),m_timeTicks(timeTicks){
  }
  
  /// Copy c'tor.  Nothing fancy, just copy all values.
  Time( const Time& other )
    :TObject(other),
     m_current(other.current()),m_previous(other.previous()),
     m_timeHack(other.timeHack()),m_timeTicks(other.timeTicks()){
  }

  /// D'tor.  Nothing special.
  virtual ~Time(){
  }
    
  /// Assignement operator.  Nothing fancy, just copy all values.
  inline Time& operator=( const Time& other ) {
    initialize(other.current(),other.previous(),
	       other.timeHack(),other.timeTicks());
    return *this;
  }
  
  /// The TimeTone right before the "active" one at event capture time
  inline const TimeTone& current() const { return m_current; } 
  
  /// The TimeTone that was "active" at event capture time
  inline const TimeTone& previous() const { return m_previous; }
  
  /// The value of the GemTime registers at event capture time
  inline const GemTime& timeHack() const { return m_timeHack; } 
  
  /// The number of 50ns ticks since last the last time hack
  inline UInt_t timeTicks() const { return m_timeTicks; }
  
  /// set everything at once
  inline void initialize(const TimeTone& current, const TimeTone& previous,
			 const GemTime& timeHack, UInt_t timeTicks) {
    m_current = current;
    m_previous = previous;      
    m_timeHack = timeHack; 
    m_timeTicks = timeTicks;
  }
  
  // set the individual data members
  inline void setCurrent( const TimeTone& val) { m_current = val; }; 
  inline void setPrevious( const TimeTone& val) { m_previous = val; };
  inline void setGemTime( const GemTime& val) { m_timeHack = val; }; 
  inline void setTimeTicks( UInt_t val) { m_timeTicks = val; };

  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_current.Clear("");
    m_previous.Clear("");
    m_timeHack.Clear("");      
    m_timeTicks = LSF_INVALID_UINT;
  }

  /// ROOT print function
  void Print(Option_t* /* option="" */) const;

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const Time &, const std::string & name = "" ) const ; // for tests
  
private:
  
  /// The TimeTone that was "active" at event capture time
  TimeTone m_current;  
  /// The TimeTone right before the "active" one at event capture time
  TimeTone m_previous;
  /// The value of the GemTime registers at event capture time
  GemTime m_timeHack;

  /// The number of 50ns ticks since last the last time hack
  UInt_t m_timeTicks;
  
  ClassDef(Time,1)   // Information about various time counters at event capture
  
};

#endif    // ROOT_TIME_H
