#ifndef LSFDATA_CONFIGURATION
#define LSFDATA_CONFIGURATION 1

#include <iostream>
#include <iomanip>
#include <string>

#include "enums/Lsf.h"

/** @class LsfEvent::Configuration
* @brief FIXME
*
* $Header$
*/

namespace lsfData {

  // forward declare all the sub-classes
  class LpaConfiguration;
  class LciConfiguration;
  class LciAcdConfiguration;
  class LciCalConfiguration;
  class LciTkrConfiguration;

  // Base configuration class, 
  // Just enough to be able to get the right type of configuration
  class Configuration {
         
  public:
    
    Configuration() {;}    
    
    virtual ~Configuration() {
    }

    virtual Configuration* clone() const {
      return 0;
    }

    virtual void clear() {
    }

    virtual void print(const std::string &str="") const {printf("%s\n",str.c_str()); }

    friend std::ostream& operator<< ( std::ostream& s, const Configuration& /*o*/) {        return s;
    }
    
    // these cast to sub-classes
    virtual const LpaConfiguration* castToLpaConfig() const { return 0; };
    virtual const LciConfiguration* castToLciConfig() const { return 0; };    
    virtual const LciAcdConfiguration* castToLciAcdConfig() const { return 0; };
    virtual const LciCalConfiguration* castToLciCalConfig() const { return 0; };
    virtual const LciTkrConfiguration* castToLciTkrConfig() const { return 0; };

    // override this to specify the Type of configuration
    virtual enums::Lsf::RunType type() const { return enums::Lsf::NoRunType; }    
  };


  // Configuration for LPA (aka particle data) running
  class LpaConfiguration : public Configuration {

  public:
    
    LpaConfiguration() 
      :Configuration(),m_hardwareKey(0),m_softwareKey(0){
    }
    
    LpaConfiguration(unsigned int hardwareKey, unsigned int softwareKey) 
      :Configuration(),m_hardwareKey(hardwareKey),m_softwareKey(softwareKey){
    }
    
    LpaConfiguration(const LpaConfiguration& other) 
      :Configuration(),m_hardwareKey(other.hardwareKey()),m_softwareKey(other.softwareKey()){
    }

    virtual ~LpaConfiguration() {
    }

    /// assignment operator
    LpaConfiguration& operator=(const LpaConfiguration& other) {
      m_hardwareKey = other.hardwareKey();
      m_softwareKey = other.softwareKey();
      return *this;
    }

    virtual Configuration* clone() const {
      return new LpaConfiguration(*this);
    }

    virtual void clear() {
        m_hardwareKey = 0;
        m_softwareKey = 0;
        Configuration::clear();
    }

    virtual void print(const std::string &str="") const {
        printf("%s softwareKey = 0x%08x\n", str.c_str(), m_softwareKey);
        printf("%s hardwareKey = 0x%08x\n", str.c_str(), m_hardwareKey);
    }

   /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const LpaConfiguration& obj )    {
      return obj.fillStream(s);
    }

    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const
 {
      return s << " softwareKey = 0x" << std::hex << std::setfill('0') 
               << m_softwareKey << "\n"
               << " hardwareKey = 0x" << m_hardwareKey << std::endl;
    }

    /// This is a poor-man's dynamic cast
    virtual const LpaConfiguration* castToLpaConfig() const { return this; };
    
    /// What type of configuration is this?
    virtual enums::Lsf::RunType type() const { return enums::Lsf::LPA; }    
    
    /// the hardware key is the LAT-C configuration master key
    inline unsigned int hardwareKey() const { return m_hardwareKey; }

    /// the software key is the key of the FSW boot file
    inline unsigned int softwareKey() const { return m_softwareKey; }

    /// set everything at once
    inline void set(unsigned int hardwareKey, unsigned int softwareKey) {
      m_hardwareKey = hardwareKey;
      m_softwareKey = softwareKey;
    }
    
    //  Set individual values
    void setHardwareKey( unsigned int value ) { m_hardwareKey = value; }
    void setSoftwareKey( unsigned int value ) { m_softwareKey = value; }  
    
    
  private:
    
    /// 
    unsigned int m_hardwareKey;
    unsigned int m_softwareKey;
    
  };


  // configuration for Lci (aka charge injection) running
  class LciConfiguration : public Configuration {  

  public:
    
    // moved to the enums package..now called by enums::Lsf::Masks
    //enum { AUTORANGE = 1, ZERO_SUPPRESSION = 2 } MASKS;

  public:
    
    LciConfiguration() 
      :m_softwareKey(0),m_writeCfg(0),m_readCfg(0),
       m_period(0),m_flags(0){
    }
     
    LciConfiguration(unsigned int softwareKey, unsigned int writeCfg, 
                     unsigned int readCfg, 
		     unsigned int period, unsigned char flags) 
      :Configuration(),
       m_softwareKey(softwareKey),m_writeCfg(writeCfg),m_readCfg(readCfg),
       m_period(period),m_flags(flags){
    }
    
    LciConfiguration(const LciConfiguration& other) 
      :Configuration(),
       m_softwareKey(other.softwareKey()),m_writeCfg(other.writeCfg()),
       m_readCfg(other.readCfg()),
       m_period(other.period()),m_flags(other.flags()){
    }

    virtual ~LciConfiguration(){
    }

    virtual void clear() {
        m_softwareKey = 0;
        m_writeCfg = 0;
        m_readCfg = 0;
        m_period = 0;
        m_flags = 0;
    }

    /// What type of configuration is this?
    virtual enums::Lsf::RunType type() const { return enums::Lsf::NoRunType; }    

    /// This is a poor-man's dynamic cast
    virtual const LciConfiguration* castToLciConfig() const { return this; };    
    /// Returns true if autorange is enabled
    inline bool autoRange() const { return (m_flags & enums::Lsf::AUTORANGE); }

    /// Returns true is zeroSuppression is on
    inline bool zeroSuppression() const { 
        return (m_flags & enums::Lsf::ZERO_SUPPRESSION);}

    /// Returns true if CALSTROBE is on
    inline bool strobe() const { return (m_flags & enums::Lsf::STROBE); }

    /// Key which uniquely identifies the LCI software configuration 
    inline unsigned int softwareKey() const { return m_softwareKey; }
    
    /// LATC file key which uniquely defines the hardward configuration 
    /// initially applied for the calibration
    inline unsigned int writeCfg() const { return m_writeCfg; }

    /// LATC file key which uniquely defines the map of LAT hardward nodes to
    /// ignore when reading the config initially applied for the calibration 
    inline unsigned int readCfg() const { return m_readCfg; }

    /// the periodic rate at which the charge injection was done (in 50ns ticks)
    inline unsigned int period() const { return m_period; }

    /// FIXME
    inline unsigned int flags() const { return m_flags; }
    
    virtual void print(const std::string &str="") const {
        printf("%s softwareKey = 0x%08x\n", str.c_str(), m_softwareKey);
        printf("%s writeCfg = 0x%08x\n", str.c_str(), m_writeCfg);
        printf("%s readCfg = 0x%08x\n", str.c_str(), m_readCfg);
        printf("%s period = %u\n", str.c_str(), m_period);
        printf("%s flags = %u\n", str.c_str(), m_flags);
    }

   /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, 
                                      const LciConfiguration& obj )    {
      return obj.fillStream(s);
    }

    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const
 {
      return s << " softwareKey = 0x" << std::hex << std::setfill('0') 
               << m_softwareKey << "\n"
               << " writeCfg = 0x" << m_writeCfg << "\n"
               << " readCfg = 0x" << m_readCfg << "\n"
               << " period = " << std::dec << m_period << "\n"
               << " flags = " << m_flags << std::endl;
             
    }

    /// set everything at once
    inline void set(unsigned int softwareKey, unsigned int writeCfg, 
                    unsigned int readCfg, 
		    unsigned int period, unsigned char flags){
      m_softwareKey = softwareKey;
      m_writeCfg = writeCfg;
      m_readCfg = readCfg;
      m_period = period;
      m_flags = flags;
    }

    // set the individual data members
    inline void setSoftwareKey( unsigned int value ) { m_softwareKey = value; }  
    inline void setWriteCfg( unsigned int value ) { m_writeCfg = value; }      
    inline void setReadCfg( unsigned int value ) { m_readCfg = value; }  
    inline void setPeriod( unsigned int value ) { m_period = value; }
    inline void setFlags( unsigned int value ) { m_flags = value; }
    inline void setAutoRange( bool value ) {
      if ( value ) m_flags |= enums::Lsf::AUTORANGE;
      else m_flags &= (~enums::Lsf::AUTORANGE);
    }
    inline void setZeroSuppression( bool value ) {
      if ( value ) m_flags |= enums::Lsf::ZERO_SUPPRESSION;
      else m_flags &= (~enums::Lsf::ZERO_SUPPRESSION);
    }
    inline void setStrobe( bool value ) {
      if ( value ) m_flags |= enums::Lsf::STROBE;
      else m_flags &= (~enums::Lsf::STROBE);
    }

          
  private:
    
    unsigned int m_softwareKey;
    unsigned int m_writeCfg;
    unsigned int m_readCfg; 
    unsigned int m_period;
    unsigned int m_flags;

  };

  class Channel {
  public:
    Channel() : m_single(0), m_all(false), m_latc(false) {};

    Channel(unsigned short single, bool all, bool latc)
     :m_single(single), m_all(all), m_latc(latc) {};

    Channel(const Channel& ch) : m_single(ch.single()), m_all(ch.all()),
            m_latc(ch.latc()) {};


    ~Channel() {}
 
    void clear() {
        m_single = 0;
        m_all = false;
        m_latc = false;
    }

    void print(const std::string &str="") const {
        printf("%s single = %u\n", str.c_str(), m_single);
        printf("%s all = %d\n", str.c_str(), m_all);
        printf("%s latc = %d\n", str.c_str(), m_latc);
    }

   /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, 
                                      const Channel& obj )    {
      return obj.fillStream(s);
    }

    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const {
      return s << " single = " << m_single << std::endl 
               << " all = " << m_all << std::endl 
               << " latc = " << m_latc << std::endl;
     }
             
    /// Returns single channel number enabled in the corresponding calibration
    /// data.  If LSF_UNDEFINED is returned, a single channel was NOT enabled
    /// and either all or latc should return a legitimate value
    inline unsigned short single() const { return m_single; }

    /// True is all channels were enabled in the corresponding calibration data
    inline bool all() const { return m_all; };

    /// Returns True if the channels enabled in the corresponding calibration
    /// data was determined by an external LATC database
    inline bool latc() const { return m_latc; };

    void set(unsigned short single, bool all, bool latc) {
        m_single = single;
        m_all = all;
        m_latc = latc;
    }

    inline void setSingle(unsigned short single) { m_single = single; };
    inline void setAll(bool all) { m_all = all; }
    inline void setLatc(bool latc) { m_latc = latc; }

  private:
     unsigned short m_single;
     bool m_all, m_latc;
  };

  // configuration for Lci (aka charge injection) ACD running
  class LciAcdConfiguration : public LciConfiguration {  

  public:

     class AcdTrigger {
     public:
         AcdTrigger() : m_veto(0), m_vernier(0), m_highDiscrim(0) {};
         AcdTrigger(unsigned short veto, unsigned short vernier,
                    unsigned short hld) : m_veto(veto), m_vernier(vernier), 
                    m_highDiscrim(hld) {};
         AcdTrigger(const AcdTrigger& trig) : m_veto(trig.veto()),
                    m_vernier(trig.vetoVernier()), 
                    m_highDiscrim(trig.highDiscrim()) {};

         ~AcdTrigger() {}

         void clear() {
             m_veto = 0;
             m_vernier = 0;
             m_highDiscrim = 0;
         }

         void print(const std::string &str="") const {
             printf("%s veto = %u\n", str.c_str(), m_veto);
             printf("%s vetoVernier = %u\n", str.c_str(), m_vernier);
             printf("%s highLevelDiscrim = %u\n", str.c_str(), m_highDiscrim);
         }

         /// Output operator (ASCII)
         friend std::ostream& operator<< ( std::ostream& s, 
                                      const AcdTrigger& obj )    {
             return obj.fillStream(s);
         }

         /// Fill the output stream (ASCII)
         std::ostream& fillStream( std::ostream& s ) const {
             return s << " veto = " << m_veto << std::endl 
                      << " vetoVernier = " << m_vernier << std::endl 
                      << " highLevelDiscrim = " << m_highDiscrim << std::endl;
         }

         /// Use in conjunction with the value returned from vetoVernier()
         /// Returns the discrimination threshold necessary to toggle ACD
         /// signal sent to the GEM.  In units of DAC counts
         /// Returns LSF_UNDEFINED if the value was determined from the LATC
         /// database.
         inline unsigned short veto() const { return m_veto; }

         /// Use in conjunction with the value returned from veto()
         /// Returns the discrimination threshold necessary to toggle ACD
         /// signal sent to the GEM.  In units of DAC counts
         /// Returns LSF_UNDEFINED if the value was determined by the LATC
         /// database.
         inline unsigned short vetoVernier() const { return m_vernier; };

         /// Returns High Level discrimination threshold necessary to toggle
         /// ACD CNO signal sent to GEM.  In units of DAC counts.
         /// Returns LSF_UNDEFINED if the value was determined by the LATC
         /// database.
         inline unsigned short highDiscrim() const { return m_highDiscrim; };

         void set(unsigned short veto, unsigned short vetoVernier,
                  unsigned short hld) { 
                  m_veto = veto;
                  m_vernier = vetoVernier;
                  m_highDiscrim = hld;
         }
  
         inline void setVeto(unsigned short veto) { m_veto = veto; }
         inline void setVetoVernier(unsigned short vetoVernier) { 
              m_vernier = vetoVernier; }
         inline void setHighDiscrim(unsigned short hld) { m_highDiscrim = hld; }

     private:
         unsigned short m_veto;
         unsigned short m_vernier;
         unsigned short m_highDiscrim;

     };

  
    LciAcdConfiguration() 
      :LciConfiguration(),
      m_injected(0),m_threshold(0),m_biasDac(0), m_holdDelay(0),
      m_hitmapDelay(0), m_range(0) {
      m_trigger.clear();
      m_channel.clear();
    }
     
    LciAcdConfiguration(unsigned short injected, unsigned short threshold,
                     unsigned short biasDac, unsigned short holdDelay, 
                     unsigned short hitmapDelay, unsigned short range,
                     AcdTrigger trigger, Channel ch ) 
      :LciConfiguration(),
       m_injected(injected),m_threshold(threshold),m_biasDac(biasDac),
       m_holdDelay(holdDelay),m_hitmapDelay(hitmapDelay), m_range(range),
       m_trigger(trigger),m_channel(ch){
    }
    
    LciAcdConfiguration(const LciAcdConfiguration& other) 
      :LciConfiguration(other.softwareKey(), other.writeCfg(),other.readCfg(),
       other.period(), other.flags()),
       m_injected(other.injected()),m_threshold(other.threshold()),
       m_biasDac(other.biasDac()), m_holdDelay(other.holdDelay()),
       m_hitmapDelay(other.hitmapDelay()), m_range(other.range()),
       m_trigger(other.trigger()),m_channel(other.channel()){
    }

    virtual ~LciAcdConfiguration(){
    }

    virtual void clear() {
        m_injected = 0;
        m_threshold = 0;
        m_biasDac = 0;
        m_holdDelay = 0;
        m_hitmapDelay = 0;
        m_range = 0;
        m_trigger.clear();
        m_channel.clear();
     }

    virtual Configuration* clone() const { return new LciAcdConfiguration( *this ); };


    /// What type of configuration is this?
    virtual enums::Lsf::RunType type() const { return enums::Lsf::AcdLCI; }    

    /// This is a poor-man's dynamic cast
    virtual const LciAcdConfiguration* castToLciAcdConfig() const { return this; };    

     virtual void print(const std::string &str="") const {
         printf("%s injected = %u\n", str.c_str(), m_injected);
         printf("%s threshold = %u\n", str.c_str(), m_threshold);
         printf("%s biasDac = %u\n", str.c_str(), m_biasDac);
         printf("%s holdDelay = %u\n", str.c_str(), m_holdDelay);
         printf("%s hitmapDelay = %u\n", str.c_str(), m_hitmapDelay);
         printf("%s range = %u\n", str.c_str(), m_range);
         m_trigger.print(str);
         m_channel.print(str);
     }

     /// Output operator (ASCII)
     friend std::ostream& operator<< ( std::ostream& s, 
                                  const LciAcdConfiguration& obj )    {
         return obj.fillStream(s);
     }

     /// Fill the output stream (ASCII)
     std::ostream& fillStream( std::ostream& s ) const {
         LciConfiguration::fillStream(s);
         return s << " injected = " << m_injected << std::endl 
                  << " threshold = " << m_threshold << std::endl 
                  << " biasDac = " << m_biasDac << std::endl
                  << " holdDelay = " << m_holdDelay << std::endl
                  << " hitmapDelay = " << m_hitmapDelay << std::endl
                  << " range = " << m_range << std::endl
                  << m_trigger
                  << m_channel << std::endl;
     }

    /// Returns the amount of charge which was injected
    /// Returns LSF_UNDEFINED if the amount was determined from the 
    /// LATC database
    inline unsigned short injected() const { return m_injected; }

    /// returns the charge threshold necessary to cross to generate the
    /// specified calibration data.  Returns LSF_UNDEFINED if determined from
    /// the LATC database
    inline unsigned short threshold() const { return m_threshold; }

    inline unsigned short biasDac() const { return m_biasDac; }
    inline unsigned short holdDelay() const { return m_holdDelay; }
    inline unsigned short hitmapDelay() const { return m_hitmapDelay; }
    inline unsigned short range() const { return m_range; }

    inline const AcdTrigger& trigger() const { return m_trigger; }
    inline const Channel& channel() const { return m_channel; }

    void set(unsigned short injected, unsigned short threshold, 
             unsigned short biasDac, unsigned short holdDelay,
             unsigned short hitmapDelay, unsigned short range,
             const AcdTrigger& trigger, const Channel& channel) {
             m_injected = injected;
             m_threshold = threshold;
             m_biasDac = biasDac;
             m_holdDelay = holdDelay;
             m_hitmapDelay = hitmapDelay;
             m_range = range;
             m_trigger = trigger;
             m_channel = channel;
    }

    void setInjected(unsigned short injected) { m_injected = injected; }
    void setThreshold(unsigned short threshold) { m_threshold = threshold; }
    void setBiasDac(unsigned short biasDac) { m_biasDac = biasDac; }
    void setHoldDelay(unsigned short holdDelay) { m_holdDelay = holdDelay; }
    void setHitmapDelay(unsigned short hitmapDelay) { m_hitmapDelay = hitmapDelay; }
    void setRange(unsigned short range) { m_range = range; }
    void setTrigger(const AcdTrigger& trig) { m_trigger = trig; }
    void setChannel(const  Channel& ch) { m_channel = ch; }
    
  private:

    unsigned short m_injected;
    unsigned short m_threshold;
    unsigned short m_biasDac;
    unsigned short m_holdDelay;
    unsigned short m_hitmapDelay;
    unsigned short m_range;
    AcdTrigger m_trigger;
    Channel    m_channel;

  };


  class LciCalConfiguration : public LciConfiguration {  

  public:

     class CalTrigger {
     public:

        CalTrigger() : m_le(0), m_he(0) {};
        CalTrigger(unsigned short le, unsigned short ITE,
                   unsigned short he, unsigned short hTE) 
        : m_le(le), m_lowTrgEna(ITE), m_he(he), m_highTrgEna(hTE) {};
        CalTrigger(const CalTrigger& cal) : m_le(cal.le()), m_he(cal.he()) {}

        ~CalTrigger() { }

        void clear() { 
            m_le = 0;
            m_lowTrgEna = 0;
            m_he = 0;
            m_highTrgEna = 0;
        }

        void print(const std::string &str="") const {
         printf("%s le = %u\n", str.c_str(), m_le);
         printf("%s lowTraEna = %u\n", str.c_str(), m_lowTrgEna);
         printf("%s he = %u\n", str.c_str(), m_he);
         printf("%s highTraEna = %u\n", str.c_str(), m_highTrgEna);
        }

        /// Output operator (ASCII)
        friend std::ostream& operator<< ( std::ostream& s, 
                                  const CalTrigger& obj )    {
            return obj.fillStream(s);
        }

        /// Fill the output stream (ASCII)
        std::ostream& fillStream( std::ostream& s ) const {
            return s << " le = " << m_le << std::endl 
                     << " lowTrgEna = " << m_lowTrgEna << std::endl
                     << " he = " << m_he << std::endl
                     << " highTrgEna = " << m_highTrgEna << std::endl;
        }

        /// Returns discrimination threshold necessary to toggle CAL's Low
        /// Energy trigger signal sent to the GEM, in units of DAC counts.
        /// Return LSF_UNDEFINED if this value was determined from the LATC
        /// database.
        inline unsigned short le() const { return m_le; };

        inline unsigned short lowTrgEna() const { return m_lowTrgEna; }

        /// Returns discrimination threshold necessary to toggle CAL's High 
        /// Energy trigger signal sent to the GEM, in units of DAC counts.
        /// Return LSF_UNDEFINED if this value was determined from the LATC
        /// database.
        inline unsigned short he() const { return m_he; };

        inline unsigned short highTrgEna() const { return m_highTrgEna; };

        void set(unsigned short le, unsigned short ITE, 
                 unsigned short he, unsigned short hTE) {
            m_le = le;
            m_lowTrgEna = ITE;
            m_he = he;
            m_highTrgEna = hTE;
        }

        inline void setLe(unsigned short le) { m_le = le; }
        inline void setLowTrgEna(unsigned short ITE) { m_lowTrgEna = ITE; };
        inline void setHe(unsigned short he) { m_he = he; }
        inline void setHighTrgEna(unsigned short hTE) { m_highTrgEna = hTE; }
   
     private:
         unsigned short m_le;
         unsigned short m_lowTrgEna;
         unsigned short m_he;
         unsigned short m_highTrgEna;
     };

    LciCalConfiguration() 
      :LciConfiguration(),
      m_uld(0),m_injected(0),m_delay(0),m_firstRange(0), m_threshold(0),
      m_calibGain(0),m_highCalEna(0),m_highRngEna(0),m_highGain(0),
      m_lowCalEna(0),m_lowRngEna(0),m_lowGain(0),
      m_trigger(),m_channel(){
    }

    LciCalConfiguration(unsigned short uld, unsigned short injected,
                        unsigned short delay, unsigned short firstRange,
                        unsigned short threshold, unsigned short calibGain,
                        unsigned short highCalEna, unsigned short highRngEna,
                        unsigned short highGain, unsigned short lowCalEna,
                        unsigned short lowRngEna, unsigned short lowGain,
                        CalTrigger trigger, Channel ch ) 
      :LciConfiguration(),
       m_uld(uld),m_injected(injected),m_delay(delay),m_firstRange(firstRange),
       m_threshold(threshold),m_calibGain(calibGain),m_highCalEna(highCalEna),
       m_highRngEna(highRngEna),m_highGain(highGain),m_lowCalEna(lowCalEna),
       m_lowRngEna(lowRngEna),m_lowGain(lowGain),
       m_trigger(trigger),m_channel(ch){
    }
    
    LciCalConfiguration(const LciCalConfiguration& other) 
      :LciConfiguration(other.softwareKey(), other.writeCfg(),other.readCfg(),
       other.period(), other.flags()),
       m_uld(other.uld()),m_injected(other.injected()),
       m_delay(other.delay()),m_firstRange(other.firstRange()),
       m_threshold(other.threshold()),
       m_calibGain(other.calibGain()), m_highCalEna(other.highCalEna()),
       m_highRngEna(other.highRngEna()), m_highGain(other.highGain()),
       m_lowCalEna(other.lowCalEna()), m_lowRngEna(other.lowRngEna()),
       m_lowGain(other.lowGain()),
       m_trigger(other.trigger()),m_channel(other.channel()){
    }

    virtual ~LciCalConfiguration(){
    }

    virtual void clear() {
        m_uld = 0;
        m_injected = 0;
        m_delay = 0;
        m_firstRange = 0;
        m_threshold = 0;
        m_calibGain = 0;
        m_highCalEna = 0;
        m_highRngEna = 0;
        m_highGain = 0;
        m_lowCalEna = 0;
        m_lowRngEna = 0;
        m_lowGain = 0;
        m_trigger.clear();
        m_channel.clear();
    }

    virtual Configuration* clone() const { return new LciCalConfiguration( *this ); };

    /// What type of configuration is this?
    virtual enums::Lsf::RunType type() const { return enums::Lsf::CalLCI; }    

    /// This is a poor-man's dynamic cast
    virtual const LciCalConfiguration* castToLciCalConfig() const { return this; };    

    virtual void print(const std::string &str="") const {
         printf("%s uld = %u\n", str.c_str(), m_uld);
         printf("%s injected = %u\n", str.c_str(), m_injected);
         printf("%s delay = %u\n", str.c_str(), m_delay);
         printf("%s firstRange = %u\n", str.c_str(), m_firstRange);
         printf("%s threshold = %u\n", str.c_str(), m_threshold);
         printf("%s calibGain = %u\n", str.c_str(), m_calibGain);
         printf("%s highCalEna = %u\n", str.c_str(), m_highCalEna);
         printf("%s highRngEna = %u\n", str.c_str(), m_highRngEna);
         printf("%s highGain = %u\n", str.c_str(), m_highGain);
         printf("%s lowCalEna = %u\n", str.c_str(), m_lowCalEna);
         printf("%s lowRngEna = %u\n", str.c_str(), m_lowRngEna);
         printf("%s lowGain = %u\n", str.c_str(), m_lowGain);
	 m_trigger.print( str.c_str() );
         m_channel.print();
     }

     /// Output operator (ASCII)
     friend std::ostream& operator<< ( std::ostream& s, 
                               const LciCalConfiguration& obj )    {
         return obj.fillStream(s);
     }

     /// Fill the output stream (ASCII)
     std::ostream& fillStream( std::ostream& s ) const {
         LciConfiguration::fillStream(s);
         return s << " uld = " << m_uld << std::endl 
                  << " injected = " << m_injected << std::endl
                  << " delay = " << m_delay << std::endl
                  << " firstRange = " << m_firstRange << std::endl
                  << " threshold = " << m_threshold << std::endl
                  << " calibGain = " << m_calibGain << std::endl
                  << " highCalEna = " << m_highCalEna << std::endl
                  << " highRngEna = " << m_highRngEna << std::endl
                  << " highGain = " << m_highGain << std::endl
                  << " lowCalEna = " << m_lowCalEna << std::endl
                  << " lowRngEna = " << m_lowRngEna << std::endl
                  << " lowGain = " << m_lowGain << std::endl
	          << m_trigger
                  << m_channel << std::endl;
      }

    /// Returns the threhsold necessary to cross range boundaries in the
    /// specified calibration data, in units of DAC counts.
    /// Returns LSF_UNDEFINED if the threshold was determined from the LATC
    /// database
    inline unsigned short uld() const { return m_uld; }

    /// Returns the amount of charge injected in units of DAC counts
    /// Returns LSF_UNDEFINED if the value was determined from the LATC
    /// database
    inline unsigned short injected() const { return m_injected; }

    /// Returns the time delay between injection of charge and TACK used to
    /// readout the corresponding calibration data.  In units of LAT clock tics
    /// where one tic is nominally 50 nanoseconds
    /// Returns LSF_UNDEFINED if the value was determined from the LATC
    /// database
    inline unsigned short delay () const { return m_delay; }

    inline unsigned short firstRange() const { return m_firstRange; }

    /// Returns threshold necessary to cross in order to generate the specified
    /// calibration data, in units of DAC counts
    /// Returns LSF_UNDEFINED if the value was determined from the LATC
    /// database
    inline unsigned short threshold() const { return m_threshold;}

    inline unsigned short calibGain() const { return m_calibGain; }
    inline unsigned short highCalEna() const { return m_highCalEna; }
    inline unsigned short highRngEna() const { return m_highRngEna; }
    inline unsigned short highGain() const { return m_highGain; }
    inline unsigned short lowCalEna() const { return m_lowCalEna; }
    inline unsigned short lowRngEna() const { return m_lowRngEna; }
    inline unsigned short lowGain() const { return m_lowGain; }

    inline const CalTrigger& trigger() const { return m_trigger; }

    /// Definition of Channel for CAL
    /// Single - Channel number in channel space.  This chanel is enabled in
    ///  all layers, in all towers.  Range [0,11]
    /// All - All channels in all FEs are enabled in all layers, over all towers
    inline const Channel& channel() const { return m_channel; }

    void set(unsigned short uld, unsigned short injected, unsigned short delay,
             unsigned short firstRange, unsigned short threshold, 
             unsigned short calibGain, unsigned short highCalEna, 
             unsigned short highRngEna, unsigned short highGain,
             unsigned short lowCalEna,
             unsigned short lowRngEna, unsigned short lowGain,
             const CalTrigger& trigger, 
             const Channel& ch) {
             m_uld = uld;
             m_injected = injected;
             m_delay = delay;
             m_firstRange = firstRange;
             m_threshold = threshold;
             m_calibGain = calibGain;
             m_highCalEna = highCalEna;
             m_highRngEna = highRngEna;
             m_highGain = highGain;
             m_lowCalEna = lowCalEna;
             m_lowRngEna = lowRngEna;
             m_lowGain = lowGain;
             m_trigger = trigger;
             m_channel = ch;
    }

    void setUld(unsigned short uld) { m_uld = uld; }
    void setInjected(unsigned short injected) { m_injected = injected; }
    void setDelay(unsigned short delay) { m_delay = delay; }
    void setFirstRange(unsigned short first) { m_firstRange = first; }
    void setThreshold(unsigned short threshold) { m_threshold = threshold; } 
    void setCalibGain(unsigned short calibGain) { m_calibGain = calibGain; }
    void setHighCalEna(unsigned short highCalEna) { m_highCalEna = highCalEna; }
    void setHighRngEna(unsigned short highRngEna) { m_highRngEna = highRngEna; }
    void setHighGain(unsigned short highGain) { m_highGain = highGain; }
    void setLowCalEna(unsigned short lowCalEna) { m_lowCalEna = lowCalEna; }
    void setLowRngEna(unsigned short lowRngEna) { m_lowRngEna = lowRngEna; }
    void setLowGain(unsigned short lowGain) { m_lowGain = lowGain; }
    void setTrigger(const CalTrigger& trig) { m_trigger = trig; }
    void setChannel(const Channel& ch) { m_channel = ch; }

  private:
      unsigned short m_uld;
      unsigned short m_injected;
      unsigned short m_delay;
      unsigned short m_firstRange;
      unsigned short m_threshold;
      unsigned short m_calibGain;
      unsigned short m_highCalEna;
      unsigned short m_highRngEna;
      unsigned short m_highGain;
      unsigned short m_lowCalEna;
      unsigned short m_lowRngEna;
      unsigned short m_lowGain;
      CalTrigger m_trigger;
      Channel m_channel;

  };


  class LciTkrConfiguration : public LciConfiguration {  

  public:
    LciTkrConfiguration() 
      :LciConfiguration(),
      m_injected(0),m_delay(0),m_threshold(0),m_splitLow(0),m_splitHigh(0),
      m_channel(){
    }

    LciTkrConfiguration(unsigned short injected, unsigned short delay,
                        unsigned short threshold, unsigned short splitLow,
                        unsigned short splitHigh, Channel ch ) 
      :LciConfiguration(),
       m_injected(injected),m_delay(delay),m_threshold(threshold),
       m_splitLow(splitLow),m_splitHigh(splitHigh),m_channel(ch){
    }
    
    LciTkrConfiguration(const LciTkrConfiguration& other) 
      :LciConfiguration(other.softwareKey(), other.writeCfg(),other.readCfg(),
       other.period(), other.flags()),
       m_injected(other.injected()),m_delay(other.delay()),
       m_threshold(other.threshold()),m_splitLow(other.splitLow()),
       m_splitHigh(other.splitHigh()), m_channel(other.channel()){
    }

    virtual ~LciTkrConfiguration(){
    }

    virtual void clear() {
        m_injected = 0;
        m_delay = 0;
        m_threshold = 0;
        m_splitLow = 0;
        m_splitHigh = 0;
        m_channel.clear();
    }

    virtual Configuration* clone() const { return new LciTkrConfiguration( *this ); };

    /// What type of configuration is this?
    virtual enums::Lsf::RunType type() const { return enums::Lsf::TkrLCI; }    

    /// This is a poor-man's dynamic cast
    virtual const LciTkrConfiguration* castToLciTkrConfig() const { return this; };    

    virtual void print(const std::string &str="") const {
         printf("%s injected = %u\n", str.c_str(), m_injected);
         printf("%s delay = %u\n", str.c_str(), m_delay);
         printf("%s threshold = %u\n", str.c_str(), m_threshold);
         printf("%s splitLow = %u\n", str.c_str(), m_splitLow);
         printf("%s splitHigh = %u\n", str.c_str(), m_splitHigh);
         m_channel.print();
     }

     /// Output operator (ASCII)
     friend std::ostream& operator<< ( std::ostream& s, 
                               const LciTkrConfiguration& obj )    {
         return obj.fillStream(s);
     }

     /// Fill the output stream (ASCII)
     std::ostream& fillStream( std::ostream& s ) const {
         LciConfiguration::fillStream(s);
         return s << " injected = " << m_injected << std::endl
                  << " delay = " << m_delay << std::endl
                  << " threshold = " << m_threshold << std::endl
                  << " splitLow = " << m_splitLow << std::endl
                  << " splitHigh = " << m_splitHigh << std::endl
                  << m_channel << std::endl;
      }

      /// Returns the amount of charge injected in units of DAC counts
      /// Returns LSF_UNDEFINED if the amount of charge injected was determined
      /// from the LATC database
      inline unsigned short injected() const { return m_injected; }


      /// Returns the time delay between injection of charge and the TACK
      /// In units of LAT clock tics where one tick is nominally 50 nanoseconds
      /// Returns LSF_UNDEFINED if the amount of charge injected was determined
      /// from the LATC database
      inline unsigned short delay() const { return m_delay; }

      /// Returns the charge threshold to cross in order to generate the 
      /// specified calibration data, in units of DAC counts.
      /// Returns LSF_UNDEFINED if the amount of charge injected was determined
      /// from the LATC database
      inline unsigned short threshold() const { return m_threshold; }

      inline unsigned short splitLow() const { return m_splitLow; }
      inline unsigned short splitHigh() const { return m_splitHigh; }

      /// TKR definition of Channel
      /// Single is the channel number in layer space, this channel is enabled
      ///  in all layers, in all towers.  Range [0,1535]
      /// All - the channel number in Front-End space.  This channel is enabled
      /// in all FEs, in all layers, in all towers.  Range [0,63]
      inline const Channel& channel() const { return m_channel;}

      void set(unsigned short injected, unsigned short delay, 
               unsigned short threshold, unsigned short splitLow,
               unsigned short splitHigh, const Channel& ch) {
               m_injected = injected;
               m_delay = delay;
               m_threshold = threshold;
               m_splitLow = splitLow;
               m_splitHigh = splitHigh;
               m_channel = ch;
      }

      inline void setInjected(unsigned short injected) { m_injected = injected; }
      inline void setDelay(unsigned short delay) { m_delay = delay; }
      inline void setThreshold(unsigned short threshold) { m_threshold = threshold; }
      inline void setSplitLow(unsigned short splitLow) { m_splitLow = splitLow; }
      inline void setSplitHigh(unsigned short splitHigh) { m_splitHigh = splitHigh; }
      inline void setChannel(const Channel& ch) { m_channel = ch; }


  private:

      unsigned short m_injected;
      unsigned short m_delay;
      unsigned short m_threshold;
      unsigned short m_splitLow;
      unsigned short m_splitHigh;
      Channel m_channel;

  };
}


#endif    // LSF_CONFIGURATION
