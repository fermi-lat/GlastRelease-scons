#ifndef ROOT_CONFIGURATION_H
#define ROOT_CONFIGURATION_H 1

#include <iostream>
#include "TObject.h"

#include "enums/Lsf.h"

/** @class LsfEvent::Configuration
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
  Configuration() {
  }    

  Configuration(const Configuration& other) 
    :TObject(other){
  }
  
  virtual ~Configuration() {
  }
  
  virtual Configuration* clone() const {
    return 0;
  }
  virtual const LpaConfiguration* castToLpaConfig() const { return 0; };
  virtual const LciConfiguration* castToLciConfig() const { return 0; };    
  virtual const LciAcdConfiguration* castToLciAcdConfig() const { return 0; };
  virtual const LciCalConfiguration* castToLciCalConfig() const { return 0; };
  virtual const LciTkrConfiguration* castToLciTkrConfig() const { return 0; };
  virtual enums::Lsf::RunType runType() const { return enums::Lsf::NoRunType; }    
  
private:

  ClassDef(Configuration,1) ;
};


class LpaConfiguration : public Configuration {
  
public:
  
  LpaConfiguration() 
    :Configuration(),m_hardwareKey(0),m_softwareKey(0){
  }
  
  LpaConfiguration(UInt_t hardwareKey, UInt_t softwareKey) 
    :Configuration(),m_hardwareKey(hardwareKey),m_softwareKey(softwareKey){
  }
  
  LpaConfiguration(const LpaConfiguration& other) 
    :Configuration(other),m_hardwareKey(other.hardwareKey()),m_softwareKey(other.softwareKey()){
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
  
  /// This is a poor-man's dynamic cast
  virtual LpaConfiguration* castToLpaConfig() { return this; };
    
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
    
  /// 
  UInt_t m_hardwareKey;
  UInt_t m_softwareKey;

  ClassDef(LpaConfiguration,1) ;  
};



class LciConfiguration : public Configuration{  

public:
  
  enum { AUTORANGE = 1,
	 ZERO_SUPPRESSION = 2 } MASKS;
  
public:
  
  LciConfiguration() 
    :m_softwareKey(0),m_writeCfg(0),m_readCfg(0),
     m_period(0),m_flags(0){
  }
  
  LciConfiguration(UInt_t softwareKey, UInt_t writeCfg, UInt_t readCfg, 
		   UInt_t period, unsigned char flags) 
    :Configuration(),
     m_softwareKey(softwareKey),m_writeCfg(writeCfg),m_readCfg(readCfg),
     m_period(period),m_flags(flags){
  }
  
  LciConfiguration(const LciConfiguration& other) 
    :Configuration(other),
     m_softwareKey(other.softwareKey()),m_writeCfg(other.writeCfg()),m_readCfg(other.readCfg()),
     m_period(other.period()),m_flags(other.flags()){
  }
  
  virtual ~LciConfiguration(){
  }
  
  /// This is a poor-man's dynamic cast
  virtual const LciConfiguration* castToLciConfig() const { return this; };    
    
  /// FIXME
  inline UInt_t softwareKey() const { return m_softwareKey; }
  
  /// FIXME
  inline UInt_t writeCfg() const { return m_writeCfg; }
  
  /// FIXME
  inline UInt_t readCfg() const { return m_readCfg; }
  
  /// the periodic rate at which the charge injection was done (in 50ns ticks)
  inline UInt_t period() const { return m_period; }
  
  /// FIXME
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
    m_softwareKey = 0;
    m_writeCfg = 0;
    m_readCfg = 0;
    m_period = 0;
    m_flags = 0;
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

  ClassDef(LciConfiguration,1);
};


#endif    // ROOT_CONFIGURATION
