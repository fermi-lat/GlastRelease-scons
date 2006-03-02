#ifndef ROOT_GEMTIME_H
#define ROOT_GEMTIME_H 1

#include "Riostream.h"
#include "TObject.h"
#include "enums/Lsf.h"

/** @class GemTime
* @brief Encapsulate the time counters stored in the Gem
*
*        Paraphrasing from the Dfi User Manual (LAT-TD-07689 v0.5, sec. 3.11)
*           Each time the 1-PPS signal arrives at the GEM, the GEM sets two values:
*             ticks -> the current value of the LAT clock, which increments at 20MHz
*             hacks -> a 7 bit counter which is incremented each 1-PPS and rolls over
*
*           These quantities are sampled at two significant times:
*             At the arrival of the 1-PPS signal at any of crates -> This sample goes in the TimeTone
*             At the event capture time.  This sample goes in the Time::timeHack() field 
*
* $Header$
*/

class GemTime : public TObject {
    
public:

  /// Default c'tor.  Assigns sentinel values to all fields  
  GemTime()
    : m_hacks(LSF_INVALID_UINT), m_ticks(LSF_INVALID_UINT){}
  
  /// Standard c'tor.  Takes input values for all fields
  GemTime( UInt_t hacks, UInt_t ticks )
    : m_hacks(hacks), m_ticks(ticks){}
  
  /// Copy c'tor.  Nothing fancy, just copy all values.
  GemTime( const GemTime& other )
    : TObject(other), m_hacks(other.hacks()), m_ticks(other.ticks()){}

  /// D'tor.  Nothing special.
  virtual ~GemTime() {
  }
  
  /// Assignement operator.  Nothing fancy, just copy all values.
  inline GemTime& operator=( const GemTime& other ) {
    initialize(other.hacks(),other.ticks());
    return *this;
  }
  
  /// number of 1-PPS time hacks.  This counter wraps around at 127)
  inline UInt_t hacks() const { return m_hacks; }
  
  /// value of the LAT system clock at the last 1-PPS time hack
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
    m_hacks = LSF_INVALID_UINT;
    m_ticks = LSF_INVALID_UINT;
  }
  /// ROOT print function
  void Print(Option_t* /* option="" */) const;

   /// Output operator (ASCII)
   friend std::ostream& operator<< ( std::ostream& s, const GemTime& obj )    {
      return obj.fillStream(s);
   }

    /// Fill the output stream (ASCII)
   std::ostream& fillStream( std::ostream& s ) const;

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const GemTime &, const std::string & name = "" ) const ; // for tests
    
private:
  
  /// Number of 1-PPS time hacks, this is a simple 7-bit counter which rolls-over
  UInt_t m_hacks;

  /// Number of ticks of 20MhZ clock at last 1-PPS time hack
  UInt_t m_ticks;  

  ClassDef(GemTime,1)        // Values of clocks latched by the Gem at 1-PPS Time Hack
  
}; 

#endif    // ROOT_GEMTIME_H
