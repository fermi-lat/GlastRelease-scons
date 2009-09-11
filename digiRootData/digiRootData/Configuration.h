#ifndef ROOT_CONFIGURATION_H
#define ROOT_CONFIGURATION_H 1

#include "Riostream.h"
#include "TObject.h"

#include "enums/Lsf.h"

/** @class Configuration
* @brief FIXME
*
* $Header$
*/

class LpaConfiguration;
class LciConfiguration;
class LciAcdConfiguration;
class LciCalConfiguration;
class LciTkrConfiguration;


class Configuration : public TObject {
  
public:

  /// Default c'tor.  Assigns sentinel values to all fields
  Configuration() {
  }    

  /// Copy c'tor.  Nothing fancy, just copy all values.
  Configuration(const Configuration& other) 
    :TObject(other){
  }
  
  /// D'tor.  Nothing special.
  virtual ~Configuration() {
  }

  /// Assignment operator.  Nothing fancy, just copy all values.
  Configuration& operator=(const Configuration& /*other*/) {
      return *this;
  }
  
  /// Virtual clone function.  Should be overridden to copy and return correct subclass
  virtual Configuration* clone() const {
    return 0;
  }

  // These are just poor man's dynamic casts
  virtual const LpaConfiguration* castToLpaConfig() const { return 0; };
  virtual const LciConfiguration* castToLciConfig() const { return 0; };    
  virtual const LciAcdConfiguration* castToLciAcdConfig() const { return 0; };
  virtual const LciCalConfiguration* castToLciCalConfig() const { return 0; };
  virtual const LciTkrConfiguration* castToLciTkrConfig() const { return 0; };

  /// return the type of configuration.  Should be overridded by sub-classes
  virtual enums::Lsf::RunType runType() const { return enums::Lsf::NoRunType; }    
  
private:

  ClassDef(Configuration,1) // Base class for LAT configuration information sent down on science data stream
};


class LpaConfiguration : public Configuration {
  
public:

  /// Default c'tor.  Assigns sentinel values to all fields
  LpaConfiguration() 
    :Configuration(),m_hardwareKey(LSF_INVALID_UINT),m_softwareKey(LSF_INVALID_UINT){
  }

  /// Standard c'tor.  Takes input values for all fields
  LpaConfiguration(UInt_t hardwareKey, UInt_t softwareKey) 
    :Configuration(),m_hardwareKey(hardwareKey),m_softwareKey(softwareKey){
  }

  /// Copy c'tor.  Nothing fancy, just copy all values.
  LpaConfiguration(const LpaConfiguration& other) 
    :Configuration(other),m_hardwareKey(other.hardwareKey()),m_softwareKey(other.softwareKey()){
  }

  /// D'tor.  Nothing special.
  virtual ~LpaConfiguration() {
  }
  
  /// Assignment operator.  Nothing fancy, just copy all values.
  LpaConfiguration& operator=(const LpaConfiguration& other) {
    m_hardwareKey = other.hardwareKey();
    m_softwareKey = other.softwareKey();
    return *this;
  }
  
  /// Override of clone function.  Copies self and returns copy.
  virtual Configuration* clone() const {
    return new LpaConfiguration(*this);
  }
  
  /// This is a poor-man's dynamic cast
  virtual const LpaConfiguration* castToLpaConfig() const { return this; };
    
  /// What type of configuration is this?
  virtual enums::Lsf::RunType runType() const { return  enums::Lsf::LPA; }    

  /// the hardware key is the LAT-C configuration master key
  inline UInt_t hardwareKey() const { return m_hardwareKey; }
  
  /// the software key is the key of the FSW boot file
  inline UInt_t softwareKey() const { return m_softwareKey; }
  
  /// set everything at once
  inline void initialize(UInt_t hardwareKey, UInt_t softwareKey) {
    m_hardwareKey = hardwareKey;
    m_softwareKey = softwareKey;
  }
  
  //  Set individual values
  void setHardwareKey( UInt_t value ) { m_hardwareKey = value; }
  void setSoftwareKey( UInt_t value ) { m_hardwareKey = value; }  
    
  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_hardwareKey = 0;
    m_softwareKey = 0;
  }

  /// ROOT print function
  void Print(Option_t* /* option="" */) const;

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const LpaConfiguration&, const std::string & name = "" ) const ; // for tests

private:
    
  /// the hardware key is the LAT-C configuration master key
  UInt_t m_hardwareKey;
  
  /// the software key is the key of the FSW boot file
  UInt_t m_softwareKey;

  ClassDef(LpaConfiguration,1) // LAT configuration information for Particle Data
};



class LciConfiguration : public Configuration{  

public:
  
  enum { AUTORANGE = 1,
	 ZERO_SUPPRESSION = 2 } MASKS;
  
public:

  /// Default c'tor.  Assigns sentinel values to all fields
  LciConfiguration() 
    :m_softwareKey(LSF_INVALID_UINT),m_writeCfg(LSF_INVALID_UINT),m_readCfg(LSF_INVALID_UINT),
     m_period(LSF_INVALID_UINT),m_flags(LSF_INVALID_UINT){
  }

  /// Standard c'tor.  Takes input values for all fields
  LciConfiguration(UInt_t softwareKey, UInt_t writeCfg, UInt_t readCfg, 
		   UInt_t period, unsigned char flags) 
    :Configuration(),
     m_softwareKey(softwareKey),m_writeCfg(writeCfg),m_readCfg(readCfg),
     m_period(period),m_flags(flags){
  }

  /// Copy c'tor.  Nothing fancy, just copy all values.
  LciConfiguration(const LciConfiguration& other) 
    :Configuration(other),
     m_softwareKey(other.softwareKey()),m_writeCfg(other.writeCfg()),m_readCfg(other.readCfg()),
     m_period(other.period()),m_flags(other.flags()){
  }
  
  /// D'tor.  Nothing special.
  virtual ~LciConfiguration(){
  }
  
  /// Assignment operator.  Nothing fancy, just copy all values.
  LciConfiguration& operator=(const LciConfiguration& other) {
      m_softwareKey = other.m_softwareKey;
      m_writeCfg = other.m_writeCfg;
      m_readCfg = other.m_readCfg;
      m_period = other.m_period;
      m_flags = other.m_flags;
      return *this;
  }

  /// Override of clone function.  Copies self and returns copy.

  virtual Configuration* clone() const {
    return new LciConfiguration(*this);
  }

  /// This is a poor-man's dynamic cast
  virtual const LciConfiguration* castToLciConfig() const { return this; };    
    
  /// the software key is the key of the FSW boot file
  inline UInt_t softwareKey() const { return m_softwareKey; }
  
  ///  FIXME
  inline UInt_t writeCfg() const { return m_writeCfg; }
  
  ///  FIXME
  inline UInt_t readCfg() const { return m_readCfg; }
  
  /// the periodic rate at which the charge injection was done (in 50ns ticks)
  inline UInt_t period() const { return m_period; }
  
  ///  FIXME
  inline UInt_t flags() const { return m_flags; }
  
  
  /// set everything at once
  inline void initialize(UInt_t softwareKey, UInt_t writeCfg, UInt_t readCfg, 
		  UInt_t period, unsigned char flags){
    m_softwareKey = softwareKey;
    m_writeCfg = writeCfg;
    m_readCfg = readCfg;
    m_period = period;
    m_flags = flags;
  }
  
  // set the individual data members
  inline void setSoftwareKey( UInt_t value ) { m_softwareKey = value; }  
  inline void setWriteCfg( UInt_t value ) { m_writeCfg = value; }      
  inline void setReadCfg( UInt_t value ) { m_readCfg = value; }  
  inline void setPeriod( UInt_t value ) { m_period = value; }
  inline void setFlags( UInt_t value ) { m_flags = value; }

  
  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_softwareKey = LSF_INVALID_UINT;
    m_writeCfg = LSF_INVALID_UINT;
    m_readCfg = LSF_INVALID_UINT;
    m_period = LSF_INVALID_UINT;
    m_flags = LSF_INVALID_UINT;
  }

  /// ROOT print function
  void Print(Option_t* /* option="" */) const;

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const LciConfiguration&, const std::string & name = "" ) const ; 

private:
  
  UInt_t m_softwareKey;
  UInt_t m_writeCfg;
  UInt_t m_readCfg;
  UInt_t m_period;
  UInt_t m_flags;

  ClassDef(LciConfiguration,1)  // Base class for configurations of charge injection data
};

class Channel : public TObject {
  public:
    Channel() : TObject(), m_single(0), m_all(false), m_latc(false) {};

    Channel(UShort_t single, Bool_t all, Bool_t latc)
     :TObject(), m_single(single), m_all(all), m_latc(latc) {};

    Channel(const Channel& ch) 
     : TObject(ch), m_single(ch.single()), m_all(ch.all()), m_latc(ch.latc()) {};


    /// Assignment operator.  Nothing fancy, just copy all values.
    Channel& operator=(const Channel& other) {
        m_single = other.m_single;
        m_all = other.m_all;
        m_latc = other.m_latc;
        return *this;
    }

    virtual ~Channel() {}
 
    void Clear(Option_t* /* option="" */) {
        m_single = 0;
        m_all = false;
        m_latc = false;
    }

    void Print(Option_t* /* option="" */) const {
        std::cout << "single: " << m_single << " all: " << m_all
                  << " latc: " << m_latc << std::endl;
    }

    /// Returns single channel number enabled in the corresponding calibration
    /// data.  If LSF_UNDEFINED is returned, a single channel was NOT enabled
    /// and either all or latc should return a legitimate value
    inline UShort_t single() const { return m_single; }

    /// True is all channels were enabled in the corresponding calibration data
    inline Bool_t all() const { return m_all; };

    /// Returns True if the channels enabled in the corresponding calibration
    /// data was determined by an external LATC database
    inline Bool_t latc() const { return m_latc; };

    void initialize(UShort_t single, Bool_t all, Bool_t latc) {
        m_single = single;
        m_all = all;
        m_latc = latc;
    }

    inline void setSingle(UShort_t single) { m_single = single; };
    inline void setAll(Bool_t all) { m_all = all; }
    inline void setLatc(Bool_t latc) { m_latc = latc; }

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const Channel&, const std::string & name = "" ) const ; // for tests

  private:
     UShort_t m_single;
     Bool_t m_all, m_latc;

     ClassDef(Channel, 1)
  };

  // configuration for Lci (aka charge injection) ACD running
  class LciAcdConfiguration : public LciConfiguration {  

  public:

     class AcdTrigger : public TObject {
     public:
         AcdTrigger() : TObject(), m_veto(0), m_vernier(0), m_highDiscrim(0) {};

         AcdTrigger(UShort_t veto, UShort_t vernier, UShort_t hld) 
           : TObject(), m_veto(veto), m_vernier(vernier), m_highDiscrim(hld) {};

         AcdTrigger(const AcdTrigger& trig) 
                    : TObject(trig), m_veto(trig.veto()),
                    m_vernier(trig.vetoVernier()), 
                    m_highDiscrim(trig.highDiscrim()) {};

         virtual ~AcdTrigger() {}
    
         /// Assignment operator.  Nothing fancy, just copy all values.
         AcdTrigger& operator=(const AcdTrigger& other) {
             m_veto = other.m_veto;
             m_vernier = other.m_vernier;
             m_highDiscrim = other.m_highDiscrim;
             return *this;
         }

         void Clear(Option_t* /* option="" */) {
             m_veto = 0;
             m_vernier = 0;
             m_highDiscrim = 0;
         }

         void Print(Option_t* /* option="" */) const {
             std::cout << "veto: " << m_veto << " vetoVernier: "  
                       << m_vernier << " highLevelDiscrim: "
                       << m_highDiscrim << std::endl;
         }

         /// Use in conjunction with the value returned from vetoVernier()
         /// Returns the discrimination threshold necessary to toggle ACD
         /// signal sent to the GEM.  In units of DAC counts
         /// Returns LSF_UNDEFINED if the value was determined from the LATC
         /// database.
         inline UShort_t veto() const { return m_veto; }

         /// Use in conjunction with the value returned from veto()
         /// Returns the discrimination threshold necessary to toggle ACD
         /// signal sent to the GEM.  In units of DAC counts
         /// Returns LSF_UNDEFINED if the value was determined by the LATC
         /// database.
         inline UShort_t vetoVernier() const { return m_vernier; };

         /// Returns High Level discrimination threshold necessary to toggle
         /// ACD CNO signal sent to GEM.  In units of DAC counts.
         /// Returns LSF_UNDEFINED if the value was determined by the LATC
         /// database.
         inline UShort_t highDiscrim() const { return m_highDiscrim; };

         void initialize(UShort_t veto, UShort_t vetoVernier,
                  UShort_t hld) { 
                  m_veto = veto;
                  m_vernier = vetoVernier;
                  m_highDiscrim = hld;
         }
  
         inline void setVeto(UShort_t veto) { m_veto = veto; }
         inline void setVetoVernier(UShort_t vetoVernier) { 
              m_vernier = vetoVernier; }
         inline void setHighDiscrim(UShort_t hld) { m_highDiscrim = hld; }

         /// Create a fake for tests
         void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

         /// Compare to another in tests
         Bool_t CompareInRange( const AcdTrigger&, const std::string & name = "" ) const ; // for tests

     private:
         UShort_t m_veto;
         UShort_t m_vernier;
         UShort_t m_highDiscrim;

         ClassDef(AcdTrigger, 1)
     };

  
    LciAcdConfiguration() 
      :LciConfiguration(),
      m_injected(0),m_threshold(0),m_biasDac(0), m_holdDelay(0),
      m_hitmapDelay(0),m_range(0) {
      m_trigger.Clear("");
      m_channel.Clear("");
    }
     
    LciAcdConfiguration(UShort_t injected, UShort_t threshold,
                     UShort_t biasDac, UShort_t holdDelay,
                     UShort_t hitmapDelay, UShort_t range, 
                     AcdTrigger trigger, Channel ch ) 
      :LciConfiguration(),
       m_injected(injected),m_threshold(threshold),m_biasDac(biasDac),
       m_holdDelay(holdDelay),m_trigger(trigger),m_channel(ch),
       m_hitmapDelay(hitmapDelay),m_range(range){
    }
    
    LciAcdConfiguration(const LciAcdConfiguration& other) 
      :LciConfiguration(other.softwareKey(), other.writeCfg(),other.readCfg(),
       other.period(), other.flags()),
       m_injected(other.injected()),m_threshold(other.threshold()),
       m_biasDac(other.biasDac()), m_holdDelay(other.holdDelay()),
       m_trigger(other.trigger()),m_channel(other.channel()),
       m_hitmapDelay(other.hitmapDelay()),m_range(other.range()){
    }

    virtual ~LciAcdConfiguration(){
    }

    /// Assignment operator.  Nothing fancy, just copy all values.
    LciAcdConfiguration& operator=(const LciAcdConfiguration& other) {
        m_injected = other.m_injected;
        m_threshold = other.m_threshold;
        m_biasDac = other.m_biasDac;
        m_holdDelay = other.m_holdDelay;
        m_trigger = other.m_trigger;
        m_channel = other.m_channel;
        m_hitmapDelay = other.m_hitmapDelay;
        m_range = other.m_range;
        LciConfiguration::initialize(other.softwareKey(), other.writeCfg(),
                          other.readCfg(),other.period(), other.flags());
        return *this;
    }


    void Clear(Option_t* /* option="" */) {
        m_injected = 0;
        m_threshold = 0;
        m_biasDac = 0;
        m_holdDelay = 0;
        m_hitmapDelay = 0;
        m_range = 0;
        m_trigger.Clear("");
        m_channel.Clear("");
     }

    virtual Configuration* clone() const { return new LciAcdConfiguration( *this ); };


    /// What type of configuration is this?
    virtual enums::Lsf::RunType type() const { return enums::Lsf::AcdLCI; }    

    /// This is a poor-man's dynamic cast
    virtual const LciAcdConfiguration* castToLciAcdConfig() const { return this; };    

    void Print(Option_t* /* option="" */) const {
         LciConfiguration::Print("");
         std::cout << "injected: " << m_injected << " threshold: " 
                   << m_threshold << " biasDac: " << m_biasDac
                   << " holdDelay: " << m_holdDelay 
                   << " hitmapDelay: " << m_hitmapDelay
                   << " range: " << m_range << std::endl;
         m_trigger.Print("");
         m_channel.Print("");
     }

    /// Returns the amount of charge which was injected
    /// Returns LSF_UNDEFINED if the amount was determined from the 
    /// LATC database
    inline UShort_t injected() const { return m_injected; }

    /// returns the charge threshold necessary to cross to generate the
    /// specified calibration data.  Returns LSF_UNDEFINED if determined from
    /// the LATC database
    inline UShort_t threshold() const { return m_threshold; }

    inline UShort_t biasDac() const { return m_biasDac; }
    inline UShort_t holdDelay() const { return m_holdDelay; }
    inline UShort_t hitmapDelay() const { return m_hitmapDelay; }
    inline UShort_t range() const { return m_range; }

    inline const AcdTrigger& trigger() const { return m_trigger; }
    inline const Channel& channel() const { return m_channel; }

    void initialize(UShort_t injected, UShort_t threshold, 
             UShort_t biasDac, UShort_t holdDelay,
             UShort_t hitmapDelay, UShort_t range,
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

    void setInjected(UShort_t injected) { m_injected = injected; }
    void setThreshold(UShort_t threshold) { m_threshold = threshold; }
    void setBiasDac(UShort_t biasDac) { m_biasDac = biasDac; }
    void setHoldDelay(UShort_t holdDelay) { m_holdDelay = holdDelay; }
    void setHitmapDelay(UShort_t hitmapDelay) { m_hitmapDelay = hitmapDelay; }
    void setRange(UShort_t range) { m_range = range; }
    void setTrigger(const AcdTrigger& trig) { m_trigger = trig; }
    void setChannel(const  Channel& ch) { m_channel = ch; }
    
    /// Create a fake for tests
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

    /// Compare to another in tests
    Bool_t CompareInRange( const LciAcdConfiguration&, const std::string & name = "" ) const ; // for tests

  private:

    UShort_t m_injected;
    UShort_t m_threshold;
    UShort_t m_biasDac;
    UShort_t m_holdDelay;
    AcdTrigger m_trigger;
    Channel    m_channel;
    UShort_t m_hitmapDelay;
    UShort_t m_range;

    ClassDef(LciAcdConfiguration,2)
  };


  class LciCalConfiguration : public LciConfiguration {  

  public:

     class CalTrigger : public TObject {
     public:

        CalTrigger() : TObject(), m_le(0), m_he(0), m_lowTrgEna(0),
                       m_highTrgEna(0) {};
        CalTrigger(unsigned short le, unsigned short lowTrgEna, 
                   unsigned short he, unsigned short highTrgEna) 
        : TObject(), m_le(le), m_he(he),  m_lowTrgEna(lowTrgEna), 
          m_highTrgEna(highTrgEna) {};
        CalTrigger(const CalTrigger& cal) 
        : TObject(cal), m_le(cal.le()), m_he(cal.he()),
          m_lowTrgEna(cal.lowTrgEna()), m_highTrgEna(cal.highTrgEna()) {}

        ~CalTrigger() { }

        /// Assignment operator.  Nothing fancy, just copy all values.
        CalTrigger& operator=(const CalTrigger& other) {
            m_le = other.m_le;
            m_lowTrgEna = other.m_lowTrgEna;
            m_he = other.m_he;
            m_highTrgEna = other.m_highTrgEna;
            return *this;
        }

        void Clear(Option_t* /* option="" */) {
            m_le = 0;
            m_lowTrgEna = 0;
            m_he = 0;
            m_highTrgEna = 0;
        }

        void Print(Option_t* /* option="" */) const {
            std::cout << "le: " << m_le << " lowTrgEna: " << m_lowTrgEna
                      << " he: " << m_he << " highTrgEna: " << m_highTrgEna
                      << std::endl;
        }

        /// Returns discrimination threshold necessary to toggle CAL's Low
        /// Energy trigger signal sent to the GEM, in units of DAC counts.
        /// Return LSF_UNDEFINED if this value was determined from the LATC
        /// database.
        inline UShort_t le() const { return m_le; };

        inline UShort_t lowTrgEna() const { return m_lowTrgEna; };

        /// Returns discrimination threshold necessary to toggle CAL's High 
        /// Energy trigger signal sent to the GEM, in units of DAC counts.
        /// Return LSF_UNDEFINED if this value was determined from the LATC
        /// database.
        inline UShort_t he() const { return m_he; };

        inline UShort_t highTrgEna() const { return m_highTrgEna; }

        void initialize(UShort_t le, UShort_t lowTrgEna, 
                        UShort_t he, UShort_t highTrgEna) {
            m_le = le;
            m_lowTrgEna = lowTrgEna;
            m_he = he;
            m_highTrgEna = highTrgEna;
        }

        inline void setLe(UShort_t le) { m_le = le; }
        inline void setLowTrgEna(UShort_t lowTrgEna) { m_lowTrgEna = lowTrgEna; }
        inline void setHe(UShort_t he) { m_he = he; }
        inline void setHighTrgEna(UShort_t highTrgEna) { m_highTrgEna = highTrgEna; }

        /// Create a fake for tests
        void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

        /// Compare to another in tests
        Bool_t CompareInRange( const CalTrigger&, const std::string & name = "" ) const ; // for tests

   
     private:
         UShort_t m_le;
         UShort_t m_he;
         UShort_t m_lowTrgEna;
         UShort_t m_highTrgEna;

         ClassDef(CalTrigger, 2)
     };

    LciCalConfiguration() 
      :LciConfiguration(),
      m_uld(0),m_injected(0),m_delay(0),m_threshold(0),
      m_trigger(),m_channel(),m_firstRange(0),m_calibGain(0),m_highCalEna(0),
      m_highRngEna(0),m_highGain(0),m_lowCalEna(0),m_lowRngEna(0),m_lowGain(0){
    }

    LciCalConfiguration(UShort_t uld, UShort_t injected,
                        UShort_t delay, UShort_t firstRange, UShort_t threshold,
                        UShort_t calibGain, UShort_t highCalEna, 
                        UShort_t highRngEna, UShort_t highGain,
                        UShort_t lowCalEna, UShort_t lowRngEna, 
                        UShort_t lowGain,
                        CalTrigger trigger, Channel ch ) 
      :LciConfiguration(),
       m_uld(uld),m_injected(injected),m_delay(delay),m_threshold(threshold),
       m_trigger(trigger),m_channel(ch),m_firstRange(firstRange),
       m_calibGain(calibGain),m_highCalEna(highCalEna),m_highRngEna(highRngEna),
       m_highGain(highGain),m_lowCalEna(lowCalEna),m_lowRngEna(lowRngEna),
       m_lowGain(lowGain){
    }
    
    LciCalConfiguration(const LciCalConfiguration& other) 
      :LciConfiguration(other.softwareKey(), other.writeCfg(),other.readCfg(),
       other.period(), other.flags()),
       m_uld(other.uld()),m_injected(other.injected()),
       m_delay(other.delay()),m_threshold(other.threshold()),
       m_trigger(other.trigger()),m_channel(other.channel()),
       m_firstRange(other.firstRange()),m_calibGain(other.calibGain()),
       m_highCalEna(other.highCalEna()),m_highRngEna(other.highRngEna()),
       m_highGain(other.highGain()),m_lowCalEna(other.lowCalEna()),
       m_lowRngEna(other.lowRngEna()),m_lowGain(other.lowGain()){
    }

    virtual ~LciCalConfiguration(){
    }

    /// Assignment operator.  Nothing fancy, just copy all values.
    LciCalConfiguration& operator=(const LciCalConfiguration& other) {
        m_uld = other.m_uld;
        m_injected = other.m_injected;
        m_delay = other.m_delay;
        m_firstRange = other.m_firstRange;
        m_threshold = other.m_threshold;
        m_calibGain = other.m_calibGain;
        m_highCalEna = other.m_highCalEna;
        m_highRngEna = other.m_highRngEna;
        m_highGain = other.m_highGain;
        m_lowCalEna = other.m_lowCalEna;
        m_lowRngEna = other.m_lowRngEna;
        m_lowGain = other.m_lowGain;
        m_trigger = other.m_trigger;
        m_channel = other.m_channel;
        LciConfiguration::initialize(other.softwareKey(), other.writeCfg(),
                        other.readCfg(),other.period(), other.flags());
        return *this;
     }

    void Clear(Option_t* /* option="" */) {
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
        m_trigger.Clear("");
        m_channel.Clear("");
    }

    virtual Configuration* clone() const { return new LciCalConfiguration( *this ); };

    /// What type of configuration is this?
    virtual enums::Lsf::RunType type() const { return enums::Lsf::CalLCI; }    

    /// This is a poor-man's dynamic cast
    virtual const LciCalConfiguration* castToLciCalConfig() const { return this; };    

    void Print(Option_t* /* option="" */) const {
         LciConfiguration::Print("");
         std::cout << "uld: " << m_uld << " injected: " << m_injected
                   << " delay: " << m_delay << " firstRange: " 
                   << m_firstRange << " threshold: " << m_threshold
                   << " calibGain: " << m_calibGain << " highCalEna: "
                   << m_highCalEna << " highRngEna: " << m_highRngEna
                   << " highGain: " << m_highGain << " lowCalEna: "
                   << m_lowCalEna << " lowRngEna: " << m_lowRngEna
                   << " lowGain: " << m_lowGain
                   << std::endl;
	 m_trigger.Print("");
         m_channel.Print("");
     }

    /// Returns the threhsold necessary to cross range boundaries in the
    /// specified calibration data, in units of DAC counts.
    /// Returns LSF_UNDEFINED if the threshold was determined from the LATC
    /// database
    inline UShort_t uld() const { return m_uld; }

    /// Returns the amount of charge injected in units of DAC counts
    /// Returns LSF_UNDEFINED if the value was determined from the LATC
    /// database
    inline UShort_t injected() const { return m_injected; }

    /// Returns the time delay between injection of charge and TACK used to
    /// readout the corresponding calibration data.  In units of LAT clock tics
    /// where one tic is nominally 50 nanoseconds
    /// Returns LSF_UNDEFINED if the value was determined from the LATC
    /// database
    inline UShort_t delay () const { return m_delay; }

    inline UShort_t firstRange () const { return m_firstRange; }

    /// Returns threshold necessary to cross in order to generate the specified
    /// calibration data, in units of DAC counts
    /// Returns LSF_UNDEFINED if the value was determined from the LATC
    /// database
    inline UShort_t threshold() const { return m_threshold;}

    inline UShort_t calibGain () const { return m_calibGain; }

    inline UShort_t highCalEna() const { return m_highCalEna; }

    inline UShort_t highRngEna() const { return m_highRngEna; }

    inline UShort_t highGain() const { return m_highGain; }

    inline UShort_t lowCalEna() const { return m_lowCalEna; }

    inline UShort_t lowRngEna() const { return m_lowRngEna; }

    inline UShort_t lowGain() const { return m_lowGain; }

    inline const CalTrigger& trigger() const { return m_trigger; }

    /// Definition of Channel for CAL
    /// Single - Channel number in channel space.  This chanel is enabled in
    ///  all layers, in all towers.  Range [0,11]
    /// All - All channels in all FEs are enabled in all layers, over all towers
    inline const Channel& channel() const { return m_channel; }

    void initialize(UShort_t uld, UShort_t injected, UShort_t delay,
             UShort_t firstRange,
             UShort_t threshold, UShort_t calibGain, UShort_t highCalEna,
             UShort_t highRngEna, UShort_t highGain, UShort_t lowCalEna,
             UShort_t lowRngEna, UShort_t lowGain,
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

    void setUld(UShort_t uld) { m_uld = uld; }
    void setInjected(UShort_t injected) { m_injected = injected; }
    void setDelay(UShort_t delay) { m_delay = delay; }
    void setFirstRange(UShort_t firstRange) { m_firstRange = firstRange; }
    void setThreshold(UShort_t threshold) { m_threshold = threshold; } 
    void setCalibGain(UShort_t calibGain) { m_calibGain = calibGain; } 
    void setHighCalEna(UShort_t highCalEna) { m_highCalEna = highCalEna; }
    void setHighRngEna(UShort_t highRngEna) { m_highRngEna = highRngEna; }
    void setHighGain(UShort_t highGain) { m_highGain = highGain; }    
    void setLowCalEna(UShort_t lowCalEna) { m_lowCalEna = lowCalEna; }
    void setLowRngEna(UShort_t lowRngEna) { m_lowRngEna = lowRngEna; }
    void setLowGain(UShort_t lowGain) { m_lowGain = lowGain; }    
    void setTrigger(const CalTrigger& trig) { m_trigger = trig; }
    void setChannel(const Channel& ch) { m_channel = ch; }

    /// Create a fake for tests
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

    /// Compare to another in tests
    Bool_t CompareInRange( const LciCalConfiguration&, const std::string & name = "" ) const ; // for tests

  private:
      UShort_t m_uld;
      UShort_t m_injected;
      UShort_t m_delay;
      UShort_t m_threshold;
      CalTrigger m_trigger;
      Channel m_channel;
      UShort_t m_firstRange;
      UShort_t m_calibGain;
      UShort_t m_highCalEna;
      UShort_t m_highRngEna;
      UShort_t m_highGain;
      UShort_t m_lowCalEna;
      UShort_t m_lowRngEna;
      UShort_t m_lowGain;

    ClassDef(LciCalConfiguration,2)
  };


  class LciTkrConfiguration : public LciConfiguration {  

  public:
    LciTkrConfiguration() 
      :LciConfiguration(),
      m_injected(0),m_delay(0),m_threshold(0),m_channel(),m_splitLow(0), 
      m_splitHigh(0) {
    }

    LciTkrConfiguration(UShort_t injected, UShort_t delay,
                        UShort_t threshold, UShort_t splitLow,
                        UShort_t splitHigh, Channel ch ) 
      :LciConfiguration(),
       m_injected(injected),m_delay(delay),m_threshold(threshold),
       m_channel(ch),m_splitLow(splitLow),m_splitHigh(splitHigh){
    }
    
    LciTkrConfiguration(const LciTkrConfiguration& other) 
      :LciConfiguration(other.softwareKey(), other.writeCfg(),other.readCfg(),
       other.period(), other.flags()),
       m_injected(other.injected()),m_delay(other.delay()),
       m_threshold(other.threshold()),m_channel(other.channel()),
       m_splitLow(other.splitLow()),m_splitHigh(other.splitHigh()){
    }

    virtual ~LciTkrConfiguration(){
    }

    /// Assignment operator.  Nothing fancy, just copy all values.
    LciTkrConfiguration& operator=(const LciTkrConfiguration& other) {
        m_injected = other.m_injected;
        m_delay = other.m_delay;
        m_threshold = other.m_threshold;
        m_splitLow = other.m_splitLow;
        m_splitHigh = other.m_splitHigh;
        m_channel = other.m_channel;
        LciConfiguration::initialize(other.softwareKey(), other.writeCfg(),
                          other.readCfg(),other.period(), other.flags());
        return *this;
    }

    void Clear(Option_t* /* option="" */) {
        m_injected = 0;
        m_delay = 0;
        m_threshold = 0;
        m_splitLow = 0;
        m_splitHigh = 0;
        m_channel.Clear("");
    }

    virtual Configuration* clone() const { return new LciTkrConfiguration( *this ); };

    /// What type of configuration is this?
    virtual enums::Lsf::RunType type() const { return enums::Lsf::TkrLCI; }    

    /// This is a poor-man's dynamic cast
    virtual const LciTkrConfiguration* castToLciTkrConfig() const { return this; };    

    void Print(Option_t* /* option="" */) const {
         LciConfiguration::Print("");
         std::cout << "injected: " << m_injected << " delay: " << m_delay
                   << " threshold: " << m_threshold << " splitLow: "
                   << m_splitLow << " splitHigh: " << m_splitHigh << std::endl;
         m_channel.Print("");
     }

      /// Returns the amount of charge injected in units of DAC counts
      /// Returns LSF_UNDEFINED if the amount of charge injected was determined
      /// from the LATC database
      inline UShort_t injected() const { return m_injected; }


      /// Returns the time delay between injection of charge and the TACK
      /// In units of LAT clock tics where one tick is nominally 50 nanoseconds
      /// Returns LSF_UNDEFINED if the amount of charge injected was determined
      /// from the LATC database
      inline UShort_t delay() const { return m_delay; }

      /// Returns the charge threshold to cross in order to generate the 
      /// specified calibration data, in units of DAC counts.
      /// Returns LSF_UNDEFINED if the amount of charge injected was determined
      /// from the LATC database
      inline UShort_t threshold() const { return m_threshold; }

      inline UShort_t splitLow() const { return m_splitLow; }
      inline UShort_t splitHigh() const { return m_splitHigh; }

      /// TKR definition of Channel
      /// Single is the channel number in layer space, this channel is enabled
      ///  in all layers, in all towers.  Range [0,1535]
      /// All - the channel number in Front-End space.  This channel is enabled
      /// in all FEs, in all layers, in all towers.  Range [0,63]
      inline const Channel& channel() const { return m_channel;}

      void initialize(UShort_t injected, UShort_t delay, 
               UShort_t threshold, UShort_t splitLow, UShort_t splitHigh,
               const Channel& ch) {
               m_injected = injected;
               m_delay = delay;
               m_threshold = threshold;
               m_splitLow = splitLow;
               m_splitHigh = splitHigh;
               m_channel = ch;
      }

      inline void setInjected(UShort_t injected) { m_injected = injected; }
      inline void setDelay(UShort_t delay) { m_delay = delay; }
      inline void setThreshold(UShort_t threshold) { m_threshold = threshold; }
      inline void setSplitLow(UShort_t splitLow) { m_splitLow = splitLow; }
      inline void setSplitHigh(UShort_t splitHigh) { m_splitHigh = splitHigh; }
      inline void setChannel(const Channel& ch) { m_channel = ch; }


    /// Create a fake for tests
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

    /// Compare to another in tests
    Bool_t CompareInRange( const LciTkrConfiguration&, const std::string & name = "" ) const ; // for tests

  private:

      UShort_t m_injected;
      UShort_t m_delay;
      UShort_t m_threshold;
      Channel m_channel;
      UShort_t m_splitLow;
      UShort_t m_splitHigh;

    ClassDef(LciTkrConfiguration,2)
  };


#endif    // ROOT_CONFIGURATION
