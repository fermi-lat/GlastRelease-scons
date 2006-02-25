#ifndef ROOT_METAEVENT_H
#define ROOT_METAEVENT_H 1

#include <iostream>
#include "TObject.h"

#include "LsfTime.h"
#include "LsfRunInfo.h"
#include "LsfDatagramInfo.h"
#include "LsfGemScalers.h"
#include "LsfConfiguration.h"

/** @class MetaEvent
* @brief Encapsulate information about the State of the LAT when a particular event was captured.
*
*        The MetaEvent consists of:
*            runInfo       -> information about the run (aka commanded acquisition) these data were take during
*            datagramInfo  -> information about the 64K datagram these data we sent down in.
*                          ->   information about autonomous mode changes are contained in the datagramInfo
*            scalers       -> extend versions of the counters in the GEM that tally data during a run
*                          ->   information about the elapsed time and deadtime are stored here
*            time          -> information about various time markers associated with the event are stored here
*            configuration -> information about configuration keys and charge injection parameters are stored here
*
* $Header$
*/

class MetaEvent : public TObject {
  
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
     m_type(configuration.runType()){
  }
 
  /// Copy c'tor.  Just copy all values.  
  /// Does a deep copy of the configuration and uses the configuration type to deduce the run type
  MetaEvent( const MetaEvent& other )
    :TObject(other),
     m_run(other.run()),
     m_datagram(other.datagram()),
     m_scalers(other.scalers()),
     m_time(other.time()),
     m_config(0),
     m_type(enums::Lsf::NoRunType){
    if ( other.configuration() != 0 ) {
      m_config = other.configuration()->clone();
      m_type = other.configuration()->runType();
    }
  }

  /// D'tor.  Delete the configuration, which had been deep-copied
  virtual ~MetaEvent(){
    delete m_config;
  }
  
  /// Information about the run this event is from
  inline const RunInfo& run() const { return m_run; }
  
  /// Information about the datagram this event came in
  inline const DatagramInfo& datagram() const { return m_datagram; }
  
  /// The extended context records
  inline const GemScalers& scalers() const { return m_scalers; }
  
  /// Information about the time markers associated with this event
  inline const Time& time() const { return m_time; } 
  
  /// Information about the configuration keys associated with this event
  inline const Configuration* configuration() const { return m_config; }

  /// Which type of run was this, particle data or charge injection 
  inline enums::Lsf::RunType runType() const { return m_type; }
  
  /// set everything at once
  inline void initialize(const RunInfo& run, const DatagramInfo& datagram, 
			 const GemScalers& scalers,
			 const Time& time,
			 const Configuration& configuration) {
    m_run = run;
    m_datagram = datagram;
    m_scalers = scalers;
    m_time = time;
    if (m_config) delete m_config;
    m_config = configuration.clone();
    m_type = configuration.runType();
  }
  
  // set the individual data members
  inline void setRun( const RunInfo& val) { m_run = val; }
  inline void setDatagram( const DatagramInfo& val) { m_datagram = val; }
  inline void setScalers( const GemScalers& val) { m_scalers = val; }
  inline void setTime( const Time& val) { m_time = val; }
  inline void setConfiguration( const Configuration& configuration ) {
    if (m_config) delete m_config;
    m_config = configuration.clone();
    m_type = configuration.runType();
  }  

  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_run.Clear("");
    m_datagram.Clear("");
    m_scalers.Clear("");
    m_time.Clear("");
    if (m_config) delete m_config;
    m_config = 0;
    m_type = enums::Lsf::NoRunType;
  }

  /// ROOT print function
  void Print(Option_t* /* option="" */) const;
  
  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 
  
  /// Compare to another in tests
  Bool_t CompareInRange( const  MetaEvent&, const std::string & name = "" ) const ; 

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
  
  ClassDef(MetaEvent,1) // information about the State of the LAT when a particular event was captured

};

#endif    // ROOT_METAEVENT
