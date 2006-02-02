#ifndef ROOT_GEMSCALARS_H
#define ROOT_GEMSCALARS_H 1

#include <iostream>
#include "TObject.h"

/** @class GemScalars
* @brief encapsulate the time tone markers that come down with the 
*
* $Header$
*/

class GemScalars : public TObject {
  
public:
  
  GemScalars( ULong64_t elapsed, ULong64_t livetime,
	      ULong64_t prescaled, ULong64_t discarded,
	      ULong64_t sequence, ULong64_t deadzone)
    :m_elapsed(elapsed),m_livetime(livetime),
     m_prescaled(prescaled),m_discarded(discarded),
     m_sequence(sequence),m_deadzone(deadzone){
  }

  GemScalars()
    :m_elapsed(0),m_livetime(0),
     m_prescaled(0),m_discarded(0),
     m_sequence(0),m_deadzone(0){
  }
  
  GemScalars(const GemScalars& other)
    :TObject(other),
     m_elapsed(other.elapsed()),m_livetime(other.livetime()),
     m_prescaled(other.prescaled()),m_discarded(other.discarded()),
     m_sequence(other.sequence()),m_deadzone(other.deadzone()){
  }
  
  virtual ~GemScalars() {}
  
  /// Assignement operator
  inline GemScalars& operator=( const GemScalars& other ) {
    initialize(other.elapsed(),other.livetime(),
	       other.prescaled(),other.discarded(),
	       other.sequence(),other.deadzone());
    return *this;
  } 
  
  /// elapsed total time since the start of the run (in 50ns ticks)
  inline ULong64_t elapsed()   const { return m_elapsed; }; 
  
  /// elapsed livetime since the start of the run (in 50ns ticks)
  inline ULong64_t livetime()  const { return m_livetime; };
  
  /// number of events prescaled away (in hardware) in this run
  inline ULong64_t prescaled() const { return m_prescaled; };
  
  /// number of events discarded (because they arrived in deadtime) in this run
  inline ULong64_t discarded() const { return m_discarded; };
  
  /// sequence number of this event within the run
  inline ULong64_t sequence()  const { return m_sequence; };
  
  /// number of events that were lost because the arrived in the dead in this run
  inline ULong64_t deadzone()  const { return m_deadzone; };
  
  /// set everything at once
  inline void initialize(ULong64_t elapsed, ULong64_t livetime,
			 ULong64_t prescaled, ULong64_t discarded,
			 ULong64_t sequence, ULong64_t deadzone) {
    m_elapsed = elapsed;
    m_livetime = livetime;
    m_prescaled = prescaled;
    m_discarded = discarded;
    m_sequence = sequence;
    m_deadzone = deadzone;
  }
  
  // set the individual data members
  inline void setElapsed( ULong64_t val ) { m_elapsed = val; }; 
  inline void setLivetime( ULong64_t val ) { m_livetime = val; };
  inline void setPrescaled( ULong64_t val ) { m_prescaled = val; };
  inline void setDiscarded( ULong64_t val ) { m_discarded = val; };
  inline void setSequence( ULong64_t val ) { m_sequence = val; };
  inline void setDeadzone( ULong64_t val ) { m_deadzone = val; };

  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_elapsed = 0;
    m_livetime = 0;
    m_prescaled = 0;
    m_discarded = 0;
    m_sequence = 0;
    m_deadzone = 0;
  }

  /// ROOT print function
  void Print(Option_t* /* option="" */) const;

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const GemScalars&, const std::string & name = "" ) const ; // for tests
    
private:
  
  ULong64_t m_elapsed;   // 
  ULong64_t m_livetime;  //
  ULong64_t m_prescaled; //
  ULong64_t m_discarded; //
  ULong64_t m_sequence;  //
  ULong64_t m_deadzone;  //

  ClassDef(GemScalars,1) ;

};


#endif    // ROOT_GEMSCALARS_H
