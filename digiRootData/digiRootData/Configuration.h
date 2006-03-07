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
  inline void set(UInt_t softwareKey, UInt_t writeCfg, UInt_t readCfg, 
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


#endif    // ROOT_CONFIGURATION
