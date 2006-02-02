#ifndef ROOT_TIMEHACK_H
#define ROOT_TIMEHACK_H 1

#include <iostream>

#include "TObject.h"

/** @class TimeHack
* @brief encapsulate the time hack markers that come down with the event context
*
* 
*
* $Header$
*/

class TimeHack : public TObject {
    
public:
  
  TimeHack()
    : m_hacks(0), m_ticks(0){}
  
  TimeHack( UInt_t hacks, UInt_t ticks )
    : m_hacks(hacks), m_ticks(ticks){}
  
  TimeHack( const TimeHack& other )
    : TObject(other), m_hacks(other.hacks()), m_ticks(other.ticks()){}
    
  ~TimeHack() {
  }
  
  /// Assignement operator
  inline TimeHack& operator=( const TimeHack& other ) {
    initialize(other.hacks(),other.ticks());
    return *this;
  }
  
  /// number of time hacks (1-second ticks) at sample time
  inline UInt_t hacks() const { return m_hacks; }
  
  /// number of time ticks (50ns second ticks) at sample time
  inline UInt_t ticks() const { return m_ticks; }
  
  /// set everything at once
  inline void initialize(UInt_t hacks, UInt_t ticks) {
    m_hacks = hacks;
    m_ticks = ticks;
  }
  
  /// Update time ticks
  inline void setTicks( UInt_t value ) { m_ticks = value; }
  
  /// Update time hacks
  inline void setHacks( UInt_t value ) { m_hacks = value; }

  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_hacks = 0;
    m_ticks = 0;
  }
  /// ROOT print function
  void Print(Option_t* /* option="" */) const;

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const TimeHack &, const std::string & name = "" ) const ; // for tests
    
private:
  
  /// Number of 1-PPS time hacks
  UInt_t m_hacks;

  /// Number of ticks of 20MhZ clock
  UInt_t m_ticks;  

  ClassDef(TimeHack,1) 
  
}; 

#endif    // ROOT_TIMEHACK_H
