#ifndef LSFDATA_METAEVENT_H
#define LSFDATA_METAEVENT_H 1

#include <iostream>

#include "lsfData/LsfTime.h"
#include "lsfData/LsfRunInfo.h"
#include "lsfData/LsfDatagramInfo.h"
#include "lsfData/LsfGemScalers.h"
#include "lsfData/LsfConfiguration.h"

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
	       const Configuration& configuration )
      :m_run(run),m_datagram(datagram),
       m_scalers(scalers),
       m_time(time),
       m_config(configuration.clone()),
       m_type(configuration.type()){
    }

    MetaEvent()
      :m_config(0),
       m_type(enums::Lsf::NoRunType){
    }
    
    MetaEvent( const MetaEvent& other ) :
       m_run(other.run()),
       m_datagram(other.datagram()),
       m_scalers(other.scalers()),
       m_time(other.time()),
       m_config(0),
       m_type(enums::Lsf::NoRunType){
      if ( other.configuration() != 0 ) {
	m_config = other.configuration()->clone();
	m_type = other.configuration()->type();
      }
    }
    
    virtual ~MetaEvent(){
      delete m_config;
    }

    inline void clear() {
       if (m_config) {
            delete m_config;
            m_config = 0;
        }
        m_run.clear();
        m_datagram.clear();
        m_scalers.clear();
        m_time.clear();      
        m_type = enums::Lsf::NoRunType;
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

    /// set everything at once
    inline void set(const RunInfo& run, const DatagramInfo& datagram, 
		    const GemScalers& scalers,
		    const Time& time,
		    const Configuration& configuration) {
      m_run = run;
      m_datagram = datagram;
      m_scalers = scalers;
      m_time = time;
      delete m_config;
      m_config = configuration.clone();
      m_type = configuration.type();
    }

    // set the individual data members
    inline void setRun( const RunInfo& val) { m_run = val; };
    inline void setDatagram( const DatagramInfo& val) { m_datagram = val; };
    inline void setScalers( const GemScalers& val) { m_scalers = val; };
    inline void setTime( const Time& val) { m_time = val; }; 
    inline void setConfiguration( const Configuration& configuration ) {
      if (m_config) delete m_config;
      m_config = configuration.clone();
      m_type = configuration.type();
    }

    
  private:
    
    /// 
    RunInfo m_run;
    DatagramInfo m_datagram;
    GemScalers m_scalers;
    Time m_time;
    Configuration* m_config;
     
    enums::Lsf::RunType m_type;

  };

}


#endif    // LSF_METAEVENT
