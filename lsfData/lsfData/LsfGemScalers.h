#ifndef LSFDATA_GEMSCALERS_H
#define LSFDATA_GEMSCALERS_H 1

#include <iostream>
#include <iomanip>
#include <cstdio>

/** @class GemScalers
* @brief FIXME
*
* $Header$
*/


namespace lsfData {

  class GemScalers {
    
  public:
    
    GemScalers( unsigned long long int elapsed, unsigned long long int livetime,
		unsigned long long int prescaled, unsigned long long int discarded,
		unsigned long long int sequence, unsigned long long int deadzone)
      :m_elapsed(elapsed),m_livetime(livetime),
       m_prescaled(prescaled),m_discarded(discarded),
       m_sequence(sequence),m_deadzone(deadzone){
    }

    GemScalers( )
      :m_elapsed(0),m_livetime(0),
       m_prescaled(0),m_discarded(0),
       m_sequence(0),m_deadzone(0){
    }

    
    GemScalers(const GemScalers& other)
      :m_elapsed(other.elapsed()),m_livetime(other.livetime()),
       m_prescaled(other.prescaled()),m_discarded(other.discarded()),
       m_sequence(other.sequence()),m_deadzone(other.deadzone()){
    }
    
    ~GemScalers() {}

    /// Assignement operator
    inline GemScalers& operator=( const GemScalers& other ) {
      set(other.elapsed(),other.livetime(),
	  other.prescaled(),other.discarded(),
	  other.sequence(),other.deadzone());
      return *this;
    }

    void clear() {
        m_elapsed = 0;
        m_livetime = 0;
        m_prescaled = 0;
        m_discarded = 0;
        m_sequence = 0;
        m_deadzone = 0;
    }

    void print() const {
      printf(" scalers:  elapsed   = 0x%016llx = %llu \n", elapsed(), elapsed());
      printf(" scalers:  livetime  = 0x%016llx = %llu \n", livetime(), livetime());
      printf(" scalers:  prescaled = 0x%016llx = %llu \n", prescaled(), prescaled());
      printf(" scalers:  discarded = 0x%016llx = %llu \n", discarded(), discarded());
      printf(" scalers:  sequence  = 0x%016llx = %llu \n", sequence(), sequence());
      printf(" scalers:  deadzone  = 0x%016llx = %llu \n", deadzone(), deadzone());
    }

  /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const GemScalers& obj )
    {
      return obj.fillStream(s);
    }

    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const {
      s << " scalers:  elapsed   = 0x" << std::hex << std::nouppercase
        << std::setfill('0') << std::setw(8)
        << elapsed() << " = " << std::dec << elapsed() << "\n"
        << " scalers:  livetime  = 0x" << std::hex << livetime() << " = "
        << std::dec << livetime() << "\n"
        << " scalers:  prescaled = 0x" << std::hex << prescaled() << " = "
        << std::dec << prescaled() << "\n"
        << " scalers:  discarded = 0x" << std::hex << discarded() << " = "
        << std::dec << discarded() << "\n"
        << " scalers:  deadzone  = 0x" << std::hex << deadzone() << " = "
        << std::dec << deadzone() << "\n";
      return s;
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
    
    /// number of events that were lost because the arrived in the dead in this run
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

    
  private:
    
    unsigned long long int m_elapsed;   // 
    unsigned long long int m_livetime;  //
    unsigned long long int m_prescaled; //
    unsigned long long int m_discarded; //
    unsigned long long int m_sequence;  //
    unsigned long long int m_deadzone;  //
   
  };

}


#endif    // LSF_GEMSCALARS_H
