#ifndef EBF_GLTCOUNTERS_H
#define EBF_GLTCOUNTERS_H


/**
 *  @class EbfGltCounters
 *  @brief Models the various counters on the GLT that are sourced
 *         for event information.
 *
 *   This models the various counters on the GLT that are sourced for
 *   event information. The most critical of these are the counters
 *   that keep track of the PPS, event and live time. There are also
 *   other counters that count window turns, events sent by the trigger
 *   to the Event Builder and events rejected because the system was
 *   busy and events rejected because of unexpired prescalers.
 *
 *   DETAILs
 *   The 3 counters are explained here. The window turn counter counts
 *   how many times the 500 nsec (approximately) coincidence window
 *   was opened. Not all opening of this window result in a trigger,
 *   and since this opening may prevent a future opening, this counter
 *   serves as input when calculating this deadtime correction.
 *   Unfortunately, there is no reasonable way to model this within
 *   the GLEAM environment.
 *   
 *   The number of events sent by the trigger is a trivial function of
 *   the number of events triggered. This is a fairly easy number to
 *   keep track of, but it does depend on one having a firm definition
 *   of what defines a trigger.
 *
 *   The number of events rejected because of unexpired prescalars is
 *   similarly a trivial number to keep track of, but with the same
 *   condition as above, i.e. one needs a firm definition of the trigger.
 *   Currently it is assumed that all triggers have no prescale so that
 *   this counter is always 0.
 *
 *   The number of events rejected because of system busy is another
 *   counter beyond the scope of being modelled in GLEAM. It is really
 *   meant as a backup to the deadtime counter. If ratio of this counter
 *   to the total number of events trigger is equal to the deadtime
 *   as measured by the ratio of 1 - livetime_counter/elapsedtime_counter,
 *   then one can be fairly well assured that the dead time is a random
 *   phenonomon (can't spell that word can we...). If however, these
 *   2 ratios are unequal, one can be certain that there is some
 *   correlation between the events taken and those dropped because
 *   of the system being busy...
 *
 *   TIME
 *   The GLT maintains 3 registers to keep track of the time
 *
 *       <li> Elapsed tick counter
 *       <li> Live tick counter
 *       <li> PPS strobe
 *
 *   These counters all count in units of the LAT system clock frequency,
 *   nominally 20MHz. The elapsed tick counter is a 25 bit register that
 *   counts continuously. When an event is triggered, the value of
 *   this counter is captured in the event tick register.
 *
 *   The live tick counter is 24 bit counter that is gated by the
 *   SYSTEM_NOT_BUSY. To avoid rollover problems, this  counter must be
 *   read at about 2Hz. This should not be a problem in normal event
 *   taking (where the event frequency is KHz.)
 *
 *   The PPS register is not really a counter, but just a register that
 *   captures the LAT system clock counter when a PPS time hack is
 *   received. This forms the basis of a venier. To get absolute time
 *   one takes the time in the event tick register and subtracts the
 *   value in this register. One multiples the difference by the
 *   system clock frequency. This frequency is self-calibrating simply
 *   by dividing at the difference of this register on after two
 *   successive 1 PPS strobes by 1 second. The absolute time is
 *   in the GPS message. The upper 7 bits of this register holds an
 *   index which can be used to correlate the current value with a
 *   GPS message.
 *
 *
 *   MODEL OF THE LAT SYSTEM CLOCK
 *   The LAT system clock is only differentially accurate. The crystal
 *   driving the clock will drift in frequency. The rate of this drift
 *   should be slow on the time scale of the 1 PPS signal. Consequently
 *   this drift can be corrected for by measuring the frequency between
 *   two successive 1 PPS signals. To make sure that we can accurately
 *   apply this correction, a model of the clock drift is built into
 *   this class.
 *
 *   The model is that the frequency of the system clock which drives
 *   the clock counter is a function of the time. For lack of a better
 *   model, this function is taken as a sinsodial variation.
 *
 *   ticks = T0 + F_nominal * T - F_nominal * A * cos (w * T) ]
 *
 *   The (however misguided) idea is that the crystal controlling the
 *   LAT system clock runs at some nominal frequency but has some
 *   orbit position variation (due to temperature) that varies the
 *   frequency.
 *
 *   F_nominal = 20 MHz   (expected frequency of the crystal)
 *   A         = 10 ** -5 (expected variation of the crystal)
 *   w         = corresponds to an orbit
 *
 *
 *   With these parameters the 
 *   ticks = 20,000,000 * T - 200 * cos (.001164 * t)
 *
 *   The maximum deviation is thus +/- 200 counts.
 *
 */

class EbfGltCounters
{
  public:

  /**
   *
   *  @fn     unsigned int ticks (double time)
   *  @brief  Returns the low 32 bits of the number of elapsed ticks
   *          of the LAT system clock
  **/
  inline unsigned int ticks(double time) const
  {
      /*
       | Haven't been careful about the absolute phasing
       | Note that 32 bits is enough to cover the necessary range. At
       | 20MHz, a 32 bit number covers about 217 seconds. We only need
       | to cover about 128 seconds.
      */
      return m_t0 + (unsigned int)(m_frequency 
                  * (time +  m_amplitude * cos (m_iperiod * time)));
  }

  
  /**
   *   @fn     unsigned int t0 (void)
   *   @brief  Returns the value of the LAT system clock at T = 0;
   *
   *    The current model of the clock variation consists of a fixed
   *    frequency term and a sinsodial variation. This parameter sets
   *    frequency of the clock. To this is added a fixed starting time
   *    or T=0 time.  A good guess would be something around 0 -
   *    the number of ticks corresponding to less than 1 second.
  **/
  inline unsigned int to (void) const { return m_t0;  }
  


  /**
   *   @fn     unsigned int t0_set (unsigned int t0) const
   *   @brief  Sets the value of the LAT system clock at T = 0;
   *
   *    The current model of the clock variation consists of a fixed
   *    frequency term and a sinsodial variation. This parameter sets
   *    frequency of the clock. To this is added a fixed starting time
   *    or T=0 time.  A good guess would be something around 0 -
   *    the number of ticks corresponding to less than 1 second.
  **/
  inline void t0_set (unsigned int t0) {  m_t0 = t0;  return;  }
      

  /**
   *   @fn     double frequency (void)
   *   @brief  Returns the frequency of the LAT system clock
   *
   *    The current model of the clock variation consists of a fixed
   *    frequency term and a sinsodial variation. This parameter sets
   *    frequency of the clock. A good guess would be something around
   *    20MHz.
  **/
  inline double frequency (void) const {  return m_frequency; }


  /**
   *   @fn     double frequency_set (double frequency)
   *   @brief  Sets the frequency of the LAT system clock
   *
   *    The current model of the clock variation consists of a fixed
   *    frequency term and a sinsodial variation. This parameter sets
   *    frequency of the clock. A good guess would be something around
   *    20MHz.
  **/
  inline void   frequency_set (double frequency)
  {
      m_frequency = frequency;
      return;
  }
      

  /**
   *   @fn     double amplitude (void) const
   *   @brief  Returns the amplitude of the variation (as a fraction)
   *
   *    The current model of the clock variation consists of a fixed
   *    frequency term and a sinsodial variation. This parameter sets
   *    period of the amplitude of the variation. A good guess would be
   *    something around 10**-5.
  **/
  inline double amplitude (void) const {  return m_amplitude; }

  /**
   *
   *   @fn    double amplitude_set (double period_in_seconds)
   *   @brief Sets the amplitude of the variation (as a fraction)
   *
   *    The current model of the clock variation consists of a fixed
   *    frequency term and a sinsodial variation. This parameter sets
   *    period of the sinsodial variation. A good guess would be
   *    something around 10**-5.
  **/
  inline void   amplitude_set (double amplitude)
  {  
      m_amplitude = amplitude;
      return;
  }

  
  /**
   *   @fn     double period (void) const
   *   @brief  Returns the period of the variation (in seconds)
   *
   *    The current model of the clock variation consists of a fixed
   *    frequency term and a sinsodial variation. This parameter sets
   *    period of the sinsodial variation. A good guess would be
   *    something corresponding to an orbit's time
  **/
  inline double period (void) const {  return  2 * 3.14159 / m_iperiod; }

  /**
   *
   *   @fn    double period_set (double period_in_seconds)
   *   @brief Sets the period of the variation (in seconds)
   *
   *    The current model of the clock variation consists of a fixed
   *    frequency term and a sinsodial variation. This parameter sets
   *    period of the sinsodial variation. A good guess would be
   *    something corresponding to an orbit's time
  **/
  inline void   period_set (double period_in_seconds)
  {
      m_iperiod = 2 * 3.14159 / period_in_seconds;
      return;
  }


  /**
   *
   *  @fn     double declare_dead (double secs)
   *  @brief  Declares this amount of dead time to be added to the
   *          accumulated deadtime counter.
   *
   *  @param  secs  The number of seconds of deadtime to declare.
   *                This is a number likely somewhat greater than 20usecs
   *  @return The accumulated deadtime in ticks.
  **/
  inline double declare_dead (double secs)
  {
      double inc = (double)(secs * m_frequency);
      return m_dead += inc;
  }

  
  /**
   *
   *  @fn     double dead (void)
   *  @brief  Returns the accumulated deadtime in clock ticks
   *
  **/
  inline double dead (void)      {   return m_dead;     }


  
  
  /**
   *
   *  @fn     unsigned int sent (void) const
   *  @brief  Returns the number of triggers sent
   *  @return The number of sent triggers
   *
   */
  inline unsigned int sent (void) const      { return m_sent;       }


  /**
   *  @fn     unsigned int increment_sent (void) const
   *  @brief  Increments the sent counter
   *  @return The number of sent triggers (inclusive of the increment)
  **/
  inline unsigned int increment_sent (void)   { return m_sent += 1;   }

  
  /**
   *
   *  @fn     double increment_sent (int count)
   *  @brief  Increments the sent count by the specified number of counts
   *
   *  @param  count The number of counts to increment the sent counter
   *  @return       The accumulated number of sent triggers
   *
  **/
  inline double increment_sent (int count)
  {
      return m_sent += count;
  }
  


  /**
   *
   *  @fn     unsigned int prescale (void) const
   *  @brief  Returns the number of events lost to the prescaler
   *  @return The number of events lost to the prescaler
   *
  **/
  inline unsigned int prescale (void) const  { return m_prescale;       }

  /**
   *
   *  @fn     unsigned int increment_prescale (void) const
   *  @brief  Increments the number of events lost to the prescaler
   *  @return The number of events lost to the prescalar
   *          (inclusive of the increment)
   *
  **/
  inline unsigned int increment_prescale (void) { return m_prescale += 1; }

  
  /**
   *
   *  @fn     unsigned int increment_prescale (int count)
   *  @brief  Increments the prescale counter by the specified number of counts
   *
   *  @param  count The number of counts to increment the sent counter
   *  @return The number of events lost to the prescalar
   *          (inclusive of the increment)
  **/
  inline unsigned int increment_prescale (int count)
  {
      return m_prescale += count;
  }

  
  /**
   *
   *  @fn     unsigned int busy (void) const
   *  @brief  Returns the number of events lost to DAQ busy
   *  @return The number of events lost to DAQ busy
   *
  **/
  inline unsigned int busy (void) const      { return m_busy;       }

  /**
   *
   *  @fn     unsigned int increment_busy (void) const
   *  @brief  Increments the number of events lost to DAQ busy
   *  @return The number of events lost to the prescalar
   *          (inclusive of the increment)
   *
  **/
  inline unsigned int increment_busy (void)   { return m_busy += 1; }

  /*!
   *
   *  @fn     double increment_prescale (int count)
   *  @brief  Increments the prescale counter by the specified number of counts
   *
   *  @param  count The number of counts to increment the busy counter
   *  @return The number of events lost to DAQ busy
   *          (inclusive of the increment)
   *
  **/
  inline unsigned int increment_busy (int count)
  {
      return m_busy += count;
  }

  
          
  /**
   *  @fn     void initialize (unsigned int t0,
   *                           unsigned int dead,
   *                           double       frequency,
   *                           double       amplitude,                       
   *                           double       period)
   *  @brief  Initializes both the counters and clock drift model
   *  @param  t0        Value to set the LAT system clock (in ticks)
   *  @param  dead      Value to set the dead time tick counter (in ticks)
   *  @param  frequency Nominal frequency of the LAT system clock (in HZ)
   *  @param  amplitude Amplitude of the LAT system clock drift
   *                    (as a fraction ie 10**-5 = 1 part in 10**5)
   *  @param  period    The period of sinsidial variations (in seconds)
   *                    (a number like 90 minutes * 60 secs/minute would
   *                     simulate an orbital frequency correction.
   *
  **/
  inline void initialize (unsigned int t0,
                          unsigned int dead,
                          double       frequency,
                          double       amplitude,                       
                          double       period)
  {
      m_t0        = t0;
      m_dead      = dead;
      m_sent      = 0;
      m_busy      = 0;
      m_prescale  = 0;
      m_frequency = frequency;
      m_amplitude = amplitude;
      m_iperiod    = 2 * 3.14159 / period;
      return;
  }
  
  inline EbfGltCounters (void )  {   return;  }

  
  /**
   *
   *   @fn    void reset (void)
   *   @brief Resets all the counters/timers
   *
  **/
  inline void reset (void)
  {
      m_t0       = 0;
      m_dead     = 0;
      m_sent     = 0;
      m_busy     = 0;
      m_prescale = 0;
  }
  
      
  
  private:

    double  m_frequency;
    /*!< The nominal clock frequency                          */
    
    double  m_iperiod;
    /*!< Clock frequency variation period                    */

    double  m_amplitude;
    /*!< Amplitude of the variation as a fraction, eg 10**-5 */

    unsigned int m_t0;
    /*!< The clock value, in ticks, at t = 0                  */

    double m_dead;
    /*!< Accumulate the total amount of deadtime */

    unsigned int m_sent;
    /*!< Number of events sent                                */

    unsigned int m_prescale;
    /*!< Number of events not taken because of prescaling     */

    unsigned int m_busy;
    /*!< Number of events not taken because of busy           */
    
};
  
     
    
#endif    
     
