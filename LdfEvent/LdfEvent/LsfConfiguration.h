#ifndef LSF_CONFIGURATION
#define LSF_CONFIGURATION 1

#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"

#include "enums/Lsf.h"
#include "lsfDataStore/LsfConfiguration.h"

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
//  class Configuration {
         
//  public:

    /// Default c'tor.
//    Configuration() {;}    
    
    /// D'tor.  Nothing special.
//    virtual ~Configuration() {
//    }

    /// Virtual clone function.  Should be overridden to copy and return correct subclass
//    virtual Configuration* clone() const {
//      return 0;
//    }
    
    // these cast to sub-classes
//    virtual const LpaConfiguration* castToLpaConfig() const { return 0; };
//    virtual const LciConfiguration* castToLciConfig() const { return 0; };    
//    virtual const LciAcdConfiguration* castToLciAcdConfig() const { return 0; };
//    virtual const LciCalConfiguration* castToLciCalConfig() const { return 0; };
//    virtual const LciTkrConfiguration* castToLciTkrConfig() const { return 0; };

    // override this to specify the Type of configuration
//    virtual enums::Lsf::RunType type() const { return enums::Lsf::NoRunType; }    
//  };


  // Configuration for LPA (aka particle data) running
  class LpaConfiguration : public lsfDataStore::LpaConfiguration {

  public:

    /// Default c'tor.  Assigns sentinel values to all fields
    LpaConfiguration() 
      :lsfDataStore::LpaConfiguration(){
    }

    /// Standard c'tor.  Takes input values for all fields
    LpaConfiguration(unsigned int hardwareKey, unsigned int softwareKey) 
      :lsfDataStore::LpaConfiguration(hardwareKey, softwareKey){
    }
    
    /// Copy c'tor.  Nothing fancy, just copy all values.
    LpaConfiguration(const LpaConfiguration& other) 
      :lsfDataStore::LpaConfiguration(other){
    }

    LpaConfiguration(const lsfDataStore::LpaConfiguration& other) 
      :lsfDataStore::LpaConfiguration(other)
    {
    }

    //LpaConfiguration(const lsfDataStore::LpaConfiguration* other) 
    //  :lsfDataStore::Configuration(),m_hardwareKey(other->hardwareKey()),
    //   m_softwareKey(other->softwareKey()) {
   // }

    /// D'tor.  Nothing special.
    virtual ~LpaConfiguration() {
    }

    /// assignment operator.  Nothing fancy, just copy all values.
    LpaConfiguration& operator=(const LpaConfiguration& other) {
      setHardwareKey(other.hardwareKey());
      setSoftwareKey(other.softwareKey());
      return *this;
    }

    /// Override of clone function.  Copies self and returns copy.
    virtual lsfDataStore::Configuration* clone() const {
      return new LpaConfiguration(*this);
    }

    /// This is a poor-man's dynamic cast
    virtual LpaConfiguration* castToLpaConfig() { return this; };
    
    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const LpaConfiguration& obj )    {
      return s << obj.hardwareKey() << ' ' << obj.softwareKey();
    }
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const LpaConfiguration& obj )    {
      return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const                            {
      return s << "class LpaConfiguration : "
	       << hardwareKey() << ' ' << softwareKey();
    }
      
    
  private:
    
  };


  // configuration for Lci (aka charge injection) running
  class LciConfiguration : public lsfDataStore::LciConfiguration {  

  public:
    
    /// Default c'tor.  Assigns sentinel values to all fields
    LciConfiguration() 
      :lsfDataStore::LciConfiguration() {
    }
     
    /// Standard c'tor.  Takes input values for all fields
    LciConfiguration(unsigned int softwareKey, unsigned int writeCfg, unsigned int readCfg, 
		     unsigned int period, unsigned char flags) 
      :lsfDataStore::LciConfiguration(softwareKey, writeCfg, readCfg, period, flags)
    {
    }

    /// Copy c'tor.  Nothing fancy, just copy all values.
    LciConfiguration(const LciConfiguration& other) 
      :lsfDataStore::LciConfiguration(other)
    {
    }

    /// D'tor.  Nothing special.
    virtual ~LciConfiguration(){
    }

    /// This is a poor-man's dynamic cast
    virtual const LciConfiguration* castToLciConfig() const { return this; };    
    
    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const LciConfiguration& obj )    {
      return s << obj.softwareKey() << ' ' << obj.writeCfg() << ' ' << obj.readCfg()
	       << ' ' << obj.period() << ' ' << obj.flags();
    }
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const LciConfiguration& obj )    {
      return obj.fillStream(s);
    }

    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const                            {
      return s << "class LciConfiguration : "
	       << softwareKey() << ' ' << writeCfg() << ' ' << readCfg() << ' ' 
	       << period() << ' ' << flags();
    }
          
  private:
    

  };

}


#endif    // LSF_CONFIGURATION
