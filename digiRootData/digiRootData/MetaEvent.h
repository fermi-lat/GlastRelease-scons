#ifndef ROOT_METAEVENT_H
#define ROOT_METAEVENT_H 1

#include "Riostream.h"
#include "TObject.h"

#include "LsfTime.h"
#include "RunInfo.h"
#include "DatagramInfo.h"
#include "GemScalers.h"
#include "Configuration.h"

/** @class MetaEvent
* @brief Encapsulate information about the State of the LAT when a particular event was captured.
*
*  The MetaEvent consists of:
*  - RunInfo  run (aka commanded acquisition) these data were take during
*  - DatagramInfo information about the 64K datagram these data we sent down in.
*                 information about autonomous mode changes 
*  - GemScalers extend versions of the counters in the GEM that tally data 
*               during a run information about the elapsed time and deadtime 
* - LsfTime various time markers associated with the event 
* - Configuration configuration keys and charge injection parameters
*
* $Header$
*/

class MetaEvent : public TObject {
  
public:

  /// Default c'tor.  Assigns sentinel or null values to all fields
  MetaEvent()
    :m_config(0),
     m_type(enums::Lsf::NoRunType),m_lpaConfig(0),m_lciAcdConfig(0),
     m_lciCalConfig(0),m_lciTkrConfig(0){
     Clear("");
  }
  
  /// Standard c'tor.  Takes input values for all fields
  MetaEvent( const RunInfo& run, const DatagramInfo& datagram, 
	     const GemScalers& scalers,
	     const LsfTime& time,
	     const Configuration& configuration )
    :m_run(run),m_datagram(datagram),
     m_scalers(scalers),
     m_time(time),
     m_config(0),m_type(enums::Lsf::NoRunType),
     m_lpaConfig(0),m_lciAcdConfig(0),m_lciCalConfig(0),
     m_lciTkrConfig(0) {
     //m_config(configuration.clone()),
     //m_type(configuration.runType()){
  }
 
  /// Copy c'tor.  Just copy all values.  
  MetaEvent( const MetaEvent& other )
    :TObject(other),
     m_run(other.run()),
     m_datagram(other.datagram()),
     m_scalers(other.scalers()),
     m_time(other.time()),
     m_type(other.runType()) {

     if (m_config) delete m_config;
     m_config = 0;
     if (m_lpaConfig) delete m_lpaConfig;
     m_lpaConfig = 0;
     if (m_lciAcdConfig) delete m_lciAcdConfig;
     m_lciAcdConfig = 0;
     if (m_lciCalConfig) delete m_lciCalConfig;
     m_lciCalConfig = 0;
     if (m_lciTkrConfig) delete m_lciTkrConfig;
     m_lciTkrConfig = 0;

     if (other.configuration())
        setConfiguration(*(other.configuration()));
     if (other.lpaConfiguration()) 
         setLpaConfiguration(*(other.lpaConfiguration()));
     if (other.lciAcdConfiguration())
         setLciAcdConfiguration(*(other.lciAcdConfiguration()));
     if (other.lciCalConfiguration())
         setLciCalConfiguration(*(other.lciCalConfiguration()));
     if (other.lciTkrConfiguration())
         setLciTkrConfiguration(*(other.lciTkrConfiguration()));


    //if ( other.configuration() != 0 ) {
    //  setConfiguration(*(other.configuration()));
      //m_config = other.configuration()->clone();
      //m_type = other.configuration()->runType();
    //}
  }

  /// D'tor.  Delete the configuration, which had been deep-copied
  virtual ~MetaEvent(){
      Clear("");
  }

 /// define assignment operator
 MetaEvent& operator=(const MetaEvent& other) {

     m_run = other.run();
     m_datagram = other.datagram();
     m_scalers = other.scalers();
     m_time = other.time();
     m_type = other.runType();

     if (m_config) delete m_config;
     m_config = 0;
     if (m_lpaConfig) delete m_lpaConfig;
     m_lpaConfig = 0;
     if (m_lciAcdConfig) delete m_lciAcdConfig;
     m_lciAcdConfig = 0;
     if (m_lciCalConfig) delete m_lciCalConfig;
     m_lciCalConfig = 0;
     if (m_lciTkrConfig) delete m_lciTkrConfig;
     m_lciTkrConfig = 0;

     if (other.configuration())
        setConfiguration(*(other.configuration()));
     if (other.lpaConfiguration())
         setLpaConfiguration(*(other.lpaConfiguration()));
     if (other.lciAcdConfiguration())
         setLciAcdConfiguration(*(other.lciAcdConfiguration()));
     if (other.lciCalConfiguration())
         setLciCalConfiguration(*(other.lciCalConfiguration()));
     if (other.lciTkrConfiguration())
         setLciTkrConfiguration(*(other.lciTkrConfiguration()));
     return *this;
  }
  
  /// Information about the run this event is from
  inline const RunInfo& run() const { return m_run; }
  
  /// Information about the datagram this event came in
  inline const DatagramInfo& datagram() const { return m_datagram; }
  
  /// The extended context records
  inline const GemScalers& scalers() const { return m_scalers; }
  
  /// Information about the time markers associated with this event
  inline const LsfTime& time() const { return m_time; } 
  
  /// Information about the configuration keys associated with this event
  inline const Configuration* configuration() const { return m_config; }

  /// Information about the configuration keys associated with this event
  inline const LpaConfiguration* lpaConfiguration() const { return m_lpaConfig; }

  /// Information about the configuration keys associated with this event
  inline const LciAcdConfiguration* lciAcdConfiguration() const { return m_lciAcdConfig; }

  /// Information about the configuration keys associated with this event
  inline const LciCalConfiguration* lciCalConfiguration() const { return m_lciCalConfig; }

  /// Information about the configuration keys associated with this event
  inline const LciTkrConfiguration* lciTkrConfiguration() const { return m_lciTkrConfig; }

  /// Which type of run was this, particle data or charge injection 
  inline enums::Lsf::RunType runType() const { return m_type; }
  
  /// set everything at once
  inline void initialize(const RunInfo& run, const DatagramInfo& datagram, 
			 const GemScalers& scalers,
			 const LsfTime& time,
			 const Configuration& configuration) {
    Clear("");
    m_run = run;
    m_datagram = datagram;
    m_scalers = scalers;
    m_time = time;
    //setConfiguration(configuration);
    //if (m_config) delete m_config;
    //m_config = configuration.clone();
    //m_type = configuration.runType();
  }
  
  // set the individual data members
  inline void setRun( const RunInfo& val) { m_run = val; }
  inline void setDatagram( const DatagramInfo& val) { m_datagram = val; }
  inline void setScalers( const GemScalers& val) { m_scalers = val; }
  inline void setLsfTime( const LsfTime& val) { m_time = val; }
  inline void setRunType( const enums::Lsf::RunType& type) { m_type = type; }
  inline void setConfiguration( const Configuration& configuration ) {
    /*if (m_config) delete m_config;
    const LpaConfiguration* ptr(0);
    switch ( configuration.runType() ) {
        case enums::Lsf::LPA:
          ptr = configuration.castToLpaConfig(); 
    }
    if (ptr) m_config = ptr->clone();
    m_type = configuration.runType();
    */
  }  

  inline void setLpaConfiguration( const LpaConfiguration& config) {
      if (m_lpaConfig) delete m_lpaConfig;
      m_lpaConfig = new LpaConfiguration(config);
  }

  inline void setLciAcdConfiguration( const LciAcdConfiguration& config) {
      if (m_lciAcdConfig) delete m_lciAcdConfig;
      m_lciAcdConfig = new LciAcdConfiguration(config);
/*
      m_lciAcdConfig = new LciAcdConfiguration();
      m_lciAcdConfig->initialize(
         config.injected(),
         config.threshold(),
         config.biasDac(), 
         config.holdDelay(), 
         config.trigger(),
         config.channel()); 
*/
      
  }

  inline void setLciCalConfiguration( const LciCalConfiguration& config) {
      if (m_lciCalConfig) delete m_lciCalConfig;
      m_lciCalConfig = new LciCalConfiguration(config);
  }

  inline void setLciTkrConfiguration( const LciTkrConfiguration& config) {
      if (m_lciTkrConfig) delete m_lciTkrConfig;
      m_lciTkrConfig = new LciTkrConfiguration(config);
  }

  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_run.Clear("");
    m_datagram.Clear("");
    m_scalers.Clear("");
    m_time.Clear("");
    if (m_config) delete m_config;
    m_config = 0;
    if (m_lpaConfig) delete m_lpaConfig;
    m_lpaConfig = 0;
    if (m_lciAcdConfig) delete m_lciAcdConfig;
    m_lciAcdConfig = 0;
    if (m_lciCalConfig) delete m_lciCalConfig;
    m_lciCalConfig = 0;
    if (m_lciTkrConfig) delete m_lciTkrConfig;
    m_lciTkrConfig = 0;
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
  LsfTime m_time;

  /// Information about the configuration keys associated with this event
  Configuration* m_config;
  
  /// Which type of run was this, particle data or charge injection 
  enums::Lsf::RunType m_type;

  LpaConfiguration *m_lpaConfig;

  LciAcdConfiguration *m_lciAcdConfig;

  LciCalConfiguration *m_lciCalConfig;

  LciTkrConfiguration *m_lciTkrConfig;
  
  ClassDef(MetaEvent,1) // information about the State of the LAT when a particular event was captured

};

#endif    // ROOT_METAEVENT
