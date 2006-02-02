#ifndef LSF_CONFIGURATION
#define LSF_CONFIGURATION 1

#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"

#include "enums/Lsf.h"

/** @class LsfEvent::Configuration
* @brief FIXME
*
* $Header$
*/

namespace LsfEvent {

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

    /// This is a poor-man's dynamic cast
    virtual LpaConfiguration* castToLpaConfig() { return this; };
    
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
    void setSoftwareKey( unsigned int value ) { m_hardwareKey = value; }  
    
    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const LpaConfiguration& obj )    {
      return s << obj.m_hardwareKey << ' ' << obj.m_softwareKey;
    }
    /// Serialize the object for reading
    friend StreamBuffer& operator>> ( StreamBuffer& s, LpaConfiguration& obj )          {
      return s >> obj.m_hardwareKey >> obj.m_softwareKey;
    }
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const LpaConfiguration& obj )    {
      return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const                            {
      return s << "class LpaConfiguration : "
	       << m_hardwareKey << ' ' << m_softwareKey;
    }
      
    
  private:
    
    /// 
    unsigned int m_hardwareKey;
    unsigned int m_softwareKey;
    
  };


  // configuration for Lci (aka charge injection) running
  class LciConfiguration : public Configuration {  

  public:
    
    enum { AUTORANGE = 1, ZERO_SUPPRESSION = 2 } MASKS;

  public:
    
    LciConfiguration() 
      :m_softwareKey(0),m_writeCfg(0),m_readCfg(0),
       m_period(0),m_flags(0){
    }
     
    LciConfiguration(unsigned int softwareKey, unsigned int writeCfg, unsigned int readCfg, 
		     unsigned int period, unsigned char flags) 
      :Configuration(),
       m_softwareKey(softwareKey),m_writeCfg(writeCfg),m_readCfg(readCfg),
       m_period(period),m_flags(flags){
    }
    
    LciConfiguration(const LciConfiguration& other) 
      :Configuration(),
       m_softwareKey(other.softwareKey()),m_writeCfg(other.writeCfg()),m_readCfg(other.readCfg()),
       m_period(other.period()),m_flags(other.flags()){
    }

    virtual ~LciConfiguration(){
    }

    /// This is a poor-man's dynamic cast
    virtual const LciConfiguration* castToLciConfig() const { return this; };    
    
    /// FIXME
    inline unsigned int softwareKey() const { return m_softwareKey; }
    
    /// FIXME
    inline unsigned int writeCfg() const { return m_writeCfg; }

    /// FIXME
    inline unsigned int readCfg() const { return m_readCfg; }

    /// the periodic rate at which the charge injection was done (in 50ns ticks)
    inline unsigned int period() const { return m_period; }

    /// FIXME
    inline unsigned int flags() const { return m_flags; }
    
    
    /// set everything at once
    inline void set(unsigned int softwareKey, unsigned int writeCfg, unsigned int readCfg, 
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

    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const LciConfiguration& obj )    {
      return s << obj.m_softwareKey << ' ' << obj.m_writeCfg << ' ' << obj.m_readCfg
	       << ' ' << obj.m_period << ' ' << obj.m_flags;
    }
    /// Serialize the object for reading
    friend StreamBuffer& operator>> ( StreamBuffer& s, LciConfiguration& obj )          {
      return s >> obj.m_softwareKey >> obj.m_writeCfg >> obj.m_readCfg >> obj.m_period >> obj.m_flags;
    }
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const LciConfiguration& obj )    {
      return obj.fillStream(s);
    }

    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const                            {
      return s << "class LciConfiguration : "
	       << m_softwareKey << ' ' << m_writeCfg << ' ' << m_readCfg << ' ' 
	       << m_period << ' ' << m_flags;
    }
          
  private:
    
    unsigned int m_softwareKey;
    unsigned int m_writeCfg;
    unsigned int m_readCfg; 
    unsigned int m_period;
    unsigned int m_flags;

  };

}


#endif    // LSF_CONFIGURATION
