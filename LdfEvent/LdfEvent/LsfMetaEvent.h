#ifndef LSF_METAEVENT_H
#define LSF_METAEVENT_H 1

#include <iostream>
#include "GaudiKernel/StreamBuffer.h"

#include "Event/TopLevel/Definitions.h"

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/DataObject.h"

#include "LdfEvent/LsfTime.h"
#include "LdfEvent/LsfRunInfo.h"
#include "LdfEvent/LsfDatagramInfo.h"
#include "LdfEvent/LsfGemScalers.h"
#include "LdfEvent/LsfConfiguration.h"

/** @class MetaEvent
*
* $Header$
*/

//static const CLID& CLID_MetaEvent = InterfaceID("MetaEvent", 1, 0);

namespace LsfEvent {
  
  class MetaEvent : public DataObject {
    
  public:
    
    /// Default c'tor.  Assigns sentinel or null values to all fields
    MetaEvent()
      :m_config(0),
       m_type(enums::Lsf::NoRunType){
    }

    /// Standard c'tor.  Takes input values for all fields
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

    /// Copy c'tor.  Just copy all values.  
    /// Does a deep copy of the configuration and uses the configuration type to deduce the run type
    MetaEvent( const MetaEvent& other )
      :DataObject(other),
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
    
    /// D'tor.  Delete the configuration, which had been deep-copied
    virtual ~MetaEvent(){
      delete m_config;
    }

    /// Retrieve reference to class definition structure
    //virtual const CLID& clID() const   { return MetaEvent::classID(); }
    //static const CLID& classID()       { return CLID_MetaEvent; }
    
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
      delete m_config;
      m_config = configuration.clone();
      m_type = configuration.type();
    }

    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const MetaEvent& obj )    {
      s << obj.m_run;
      s << obj.m_datagram;
      s << obj.m_scalers;
      s << obj.m_time;
      s << obj.m_type; 
      if ( obj.m_config != 0) {
	// FIXME cast to correct type
	// s << *obj.m_config;
      };
      return s;
    }
    /// Serialize the object for reading
    friend StreamBuffer& operator>> ( StreamBuffer& s, MetaEvent& obj )          {
      s >> obj.m_run;
      s >> obj.m_datagram;
      s >> obj.m_scalers;
      s >> obj.m_time;
      // FIXME can't read enums
      // s >> obj.m_type;
      obj.m_config = 0;
      switch ( obj.m_type ) {
      case enums::Lsf::NoRunType: break;
      case enums::Lsf::LPA: 
	obj.m_config = new LpaConfiguration(); break;
      case enums::Lsf::AcdLCI: 
      case enums::Lsf::CalLCI:
      case enums::Lsf::TkrLCI:
	break;
      }
      if ( obj.m_config != 0 ) {
	// FIXME cast to correct type
	// s >> *obj.m_config;
      }
      return s;
    }
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const MetaEvent& obj )    {
      return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const {
      s << "class MetaEvent : "
	<< EventField( EventFormat::field12 )
	<< m_run 
	<< m_datagram 
	<< m_scalers 
	<< m_time 
	<< m_type;
      if ( m_config != 0 ) {
	// FIXME cast to correct type
	// s << *m_config 
      }
      return s;
    }
    
  private:
    
    /// Information about the run this event is from
    RunInfo m_run;
    
    /// Information about the datagram this event came in
    DatagramInfo m_datagram;
    
    /// The extended context records
    GemScalers m_scalers;

    /// Information about the time markers associated with this event  
    Time m_time;
    
    /// Information about the configuration keys associated with this event
    Configuration* m_config;    //-> 
    
    /// Which type of run was this, particle data or charge injection 
    enums::Lsf::RunType m_type;

  };

}


#endif    // LSF_METAEVENT
