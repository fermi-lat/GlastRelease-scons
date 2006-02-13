#ifndef LSF_GEMSCALERS_H
#define LSF_GEMSCALERS_H 1

#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"

#include "enums/Lsf.h"

/** @class GemScalers
* @brief FIXME
*
* $Header$
*/

//static const CLID& CLID_LdfTime = InterfaceID("GemScalers", 1, 0);

namespace LsfEvent {

  class GemScalers {
    
  public:

    /// Default c'tor.  Assigns sentinel values to all fields
    GemScalers()
      :m_elapsed(LSF_INVALID_ULONG),m_livetime(LSF_INVALID_ULONG),
       m_prescaled(LSF_INVALID_ULONG),m_discarded(LSF_INVALID_ULONG),
       m_sequence(LSF_INVALID_ULONG),m_deadzone(LSF_INVALID_ULONG){
    }
    
    /// Standard c'tor.  Takes input values for all fields  
    GemScalers( unsigned long long int elapsed, unsigned long long int livetime,
		unsigned long long int prescaled, unsigned long long int discarded,
		unsigned long long int sequence, unsigned long long int deadzone)
      :m_elapsed(elapsed),m_livetime(livetime),
       m_prescaled(prescaled),m_discarded(discarded),
       m_sequence(sequence),m_deadzone(deadzone){
    }
    
    /// Copy c'tor.  Nothing fancy, just copy all values.
    GemScalers(const GemScalers& other)
      :m_elapsed(other.elapsed()),m_livetime(other.livetime()),
       m_prescaled(other.prescaled()),m_discarded(other.discarded()),
       m_sequence(other.sequence()),m_deadzone(other.deadzone()){
    }
    
    /// D'tor.  Nothing special.
    virtual ~GemScalers() {}
    
    /// Assignement operator.  Nothing fancy, just copy all values.
    inline GemScalers& operator=( const GemScalers& other ) {
      set(other.elapsed(),other.livetime(),
	  other.prescaled(),other.discarded(),
	  other.sequence(),other.deadzone());
      return *this;
    } 
    
    /// elapsed total time since the start of the run (in 50ns ticks)
    inline unsigned long long int elapsed()   const { return m_elapsed; }; 
    
    /// elapsed livetime since the start of the run (in 50ns ticks)
    inline unsigned long long int livetime()  const { return m_livetime; };
    
    /// number of events prescaled away (in hardware) in this run
    inline unsigned long long int prescaled() const { return m_prescaled; };
    
    /// number of events discarded (because they arrived in deadtime) in this run
    inline unsigned long long int discarded() const { return m_discarded; };
    
    /// sequence number of this event within the run
    inline unsigned long long int sequence()  const { return m_sequence; };
    
    /// number of events that were lost because the arrived in the 2-tick deadzone in this run
    inline unsigned long long int deadzone()  const { return m_deadzone; };
    
    /// set everything at once
    inline void set(unsigned long long int elapsed, unsigned long long int livetime,
		    unsigned long long int prescaled, unsigned long long int discarded,
		    unsigned long long int sequence, unsigned long long int deadzone) {
      m_elapsed = elapsed;
      m_livetime = livetime;
      m_prescaled = prescaled;
      m_discarded = discarded;
      m_sequence = sequence;
      m_deadzone = deadzone;
    }
    
    // set the individual data members
    inline void setElapsed( unsigned long long int val ) { m_elapsed = val; }; 
    inline void setLivetime( unsigned long long int val ) { m_livetime = val; };
    inline void setPrescaled( unsigned long long int val ) { m_prescaled = val; };
    inline void setDiscarded( unsigned long long int val ) { m_discarded = val; };
    inline void setSequence( unsigned long long int val ) { m_sequence = val; };
    inline void setDeadzone( unsigned long long int val ) { m_deadzone = val; };
    
    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const GemScalers& obj )    {
      // FIXME StreamBuffer can't handle unsigned long long int
      //return s << obj.m_elapsed << ' ' << obj.m_livetime << ' ' << obj.m_prescaled << ' ' 
      // << obj.m_discarded << ' ' << obj.m_sequence << ' ' << obj.m_deadzone;
      return s;
    }
    /// Serialize the object for reading
    friend StreamBuffer& operator>> ( StreamBuffer& s, GemScalers& obj )          {
      // FIXME StreamBuffer can't handle unsigned long long int
      //return s >> obj.m_elapsed >> obj.m_livetime >> obj.m_prescaled 
      //       >> obj.m_discarded >> obj.m_sequence >> obj.m_deadzone;
      return s;
    }
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const GemScalers& obj )    {
      return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const                            {
      return s << "class GemScalers : "
	       << EventField( EventFormat::field12 )
	       << m_elapsed << ' ' << m_livetime << ' ' << m_prescaled << ' ' 
	       << m_discarded << ' ' << m_sequence << ' ' << m_deadzone;
    }
      
    
  private:
    
    /// elapsed total time since the start of the run (in 50ns ticks)
    unsigned long long int m_elapsed;   
    /// elapsed livetime since the start of the run (in 50ns ticks)
    unsigned long long int m_livetime; 
    /// number of events prescaled away (in hardware) in this run
    unsigned long long int m_prescaled;
    /// number of events discarded (because they arrived in deadtime) in this run
    unsigned long long int m_discarded; 
    /// sequence number of this event within the run
    unsigned long long int m_sequence;
    /// number of events that were lost because the arrived in the 2-tick deadzone in this run
    unsigned long long int m_deadzone;
   
  };

}


#endif    // LSF_GEMSCALERS_H
