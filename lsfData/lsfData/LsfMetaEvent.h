#ifndef LSFDATA_METAEVENT_H
#define LSFDATA_METAEVENT_H 1

#include <iostream>

#include "lsfData/LsfTime.h"
#include "lsfData/LsfRunInfo.h"
#include "lsfData/LsfDatagramInfo.h"
#include "lsfData/LsfGemScalers.h"
#include "lsfData/LsfConfiguration.h"
#include "lsfData/LsfKeys.h"

/** @class MetaEvent
*
* $Header$
*/

namespace lsfData {
  
  class MetaEvent {
    
  public:
    
    MetaEvent( const RunInfo& run, const DatagramInfo& datagram, 
	       const GemScalers& scalers,
	       const Time& time,
	       const Configuration& configuration,
	       const LsfKeys& keys)
      :m_run(run),m_datagram(datagram),
       m_scalers(scalers),
       m_time(time),
       m_config(configuration.clone()),
       m_type(configuration.type()),
       m_keys(keys.clone()),
       m_ktype(keys.type()){
    }

    MetaEvent()
      :m_config(0),
       m_type(enums::Lsf::NoRunType),
       m_keys(0),
       m_ktype(enums::Lsf::NoKeysType){
    }
    
    MetaEvent( const MetaEvent& other ) :
       m_run(other.run()),
       m_datagram(other.datagram()),
       m_scalers(other.scalers()),
       m_time(other.time()),
       m_config(0),
       m_type(enums::Lsf::NoRunType) {
      if ( other.configuration() != 0 ) {
	m_config = other.configuration()->clone();
	m_type = other.configuration()->type();
      }
      if ( other.keys() != 0 ) {
	m_keys  = other.keys()->clone();
	m_ktype = other.keys()->type();
      }
    }
    
    virtual ~MetaEvent(){
      delete m_config;
      delete m_keys;
    }

    inline void clear() {
       if (m_config) {
            delete m_config;
            m_config = 0;
        }
       if (m_keys) {
	 delete m_keys;
	 m_keys = 0;
       }
        m_run.clear();
        m_datagram.clear();
        m_scalers.clear();
        m_time.clear();      
        m_type = enums::Lsf::NoRunType;
	m_ktype = enums::Lsf::NoKeysType;
    }

    /// Information about the run this event is from
    inline const RunInfo& run() const { return m_run; };

    /// Information about the datagram this event came in
    inline const DatagramInfo& datagram() const { return m_datagram; }

    /// The extended context records
    inline const GemScalers& scalers() const { return m_scalers; }

    /// Information about the time markers associated with this event
    inline const Time& time() const { return m_time; } 

    /// Information about the configuration keys associated with this event
    inline const Configuration* configuration() const { return m_config; }

    /// Translated configuration file keys for this event
    inline const LsfKeys* keys() const { return m_keys; };

    /// set everything at once
    inline void set(const RunInfo& run, const DatagramInfo& datagram, 
		    const GemScalers& scalers,
		    const Time& time,
		    const Configuration& configuration,
		    const LsfKeys& keys) {
      m_run = run;
      m_datagram = datagram;
      m_scalers = scalers;
      m_time = time;
      if(m_config) delete m_config;
      m_config = configuration.clone();
      m_type = configuration.type();
      if(m_keys) delete m_keys;
      m_keys = keys.clone();
      m_ktype = keys.type();
    }

    // set the individual data members
    inline void setRun( const RunInfo& val) { m_run = val; };
    inline void setRunDataTransferId(unsigned int id) { m_run.setDataTransferId(id); };
    inline void setDatagram( const DatagramInfo& val) { m_datagram = val; };
    inline void setScalers( const GemScalers& val) { m_scalers = val; };
    inline void setTime( const Time& val) { m_time = val; }; 
    inline void setConfiguration( const Configuration& configuration ) {
      if (m_config) delete m_config;
      m_config = configuration.clone();
      m_type = configuration.type();
    }
    inline void setKeys( const LsfKeys& keys ) {
      if (m_keys) delete m_keys;
      m_keys = keys.clone();
      m_ktype = keys.type();
    }
    
  private:
    
    /// 
    RunInfo m_run;
    DatagramInfo m_datagram;
    GemScalers m_scalers;
    Time m_time;
    Configuration* m_config;
     
    enums::Lsf::RunType m_type;
    
    LsfKeys* m_keys;
    enums::Lsf::KeysType m_ktype;

  };

}


#endif    // LSF_METAEVENT
