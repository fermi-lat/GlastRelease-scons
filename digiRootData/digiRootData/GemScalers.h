#ifndef ROOT_GEMSCALERS_H
#define ROOT_GEMSCALERS_H 1

#include <iostream>
#include "TObject.h"

#include "enums/Lsf.h"

/** @class GemScalers
* @brief Encapsulate the extended counter versions held by FSW of the GemScalers 
*
*        Six counters are passed down:
*            elapsed     -> the number of 50ns ticks since the start of the run
*            livetime    -> the number of 50ns ticks that the Gem was enabled since the start of the run
*            prescaled   -> the total number of events prescaled away up to this point in the run
*            discarded   -> the total number of event discarded due to deadtime since the start of the run
*            sequence    -> the absolute sequence number of the event within the run as seen by the Gem
*            deadzone    -> the total number of event discarded b/c they arrived during the 2tick Gem deadzone
*         
* $Header$
*/

class GemScalers : public TObject {
  
public:

  /// Default c'tor.  Assigns sentinel values to all fields
  GemScalers()
    :m_elapsed(LSF_INVALID_ULONG),m_livetime(LSF_INVALID_ULONG),
     m_prescaled(LSF_INVALID_ULONG),m_discarded(LSF_INVALID_ULONG),
     m_sequence(LSF_INVALID_ULONG),m_deadzone(LSF_INVALID_ULONG){
  }
  
  /// Standard c'tor.  Takes input values for all fields  
  GemScalers( ULong64_t elapsed, ULong64_t livetime,
	      ULong64_t prescaled, ULong64_t discarded,
	      ULong64_t sequence, ULong64_t deadzone)
    :m_elapsed(elapsed),m_livetime(livetime),
     m_prescaled(prescaled),m_discarded(discarded),
    m_sequence(sequence),m_deadzone(deadzone){
  }

  /// Copy c'tor.  Nothing fancy, just copy all values.
  GemScalers(const GemScalers& other)
    :TObject(other),
     m_elapsed(other.elapsed()),m_livetime(other.livetime()),
     m_prescaled(other.prescaled()),m_discarded(other.discarded()),
     m_sequence(other.sequence()),m_deadzone(other.deadzone()){
  }

  /// D'tor.  Nothing special.
  virtual ~GemScalers() {}
  
  /// Assignement operator.  Nothing fancy, just copy all values.
  inline GemScalers& operator=( const GemScalers& other ) {
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
  
  /// number of events that were lost because the arrived in the 2-tick deadzone in this run
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
    m_elapsed = LSF_INVALID_ULONG;
    m_livetime = LSF_INVALID_ULONG;
    m_prescaled = LSF_INVALID_ULONG;
    m_discarded = LSF_INVALID_ULONG;
    m_sequence = LSF_INVALID_ULONG;
    m_deadzone = LSF_INVALID_ULONG;
  }

  /// ROOT print function
  void Print(Option_t* /* option="" */) const;

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const GemScalers&, const std::string & name = "" ) const ; // for tests
    
private:

  /// elapsed total time since the start of the run (in 50ns ticks)
  ULong64_t m_elapsed;   
  /// elapsed livetime since the start of the run (in 50ns ticks)
  ULong64_t m_livetime; 
  /// number of events prescaled away (in hardware) in this run
  ULong64_t m_prescaled;
  /// number of events discarded (because they arrived in deadtime) in this run
  ULong64_t m_discarded; 
  /// sequence number of this event within the run
  ULong64_t m_sequence;
  /// number of events that were lost because the arrived in the 2-tick deadzone in this run
  ULong64_t m_deadzone;

  ClassDef(GemScalers,1)  // Information about the FSW extended verisons of Gem counters

};


#endif    // ROOT_GEMSCALERS_H
